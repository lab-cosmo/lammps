/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_kokkos.h"

#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr int EXTRA = 1000;

/* ----------------------------------------------------------------------
   allocate and initialize array or hash table for global -> local map
   for array option:
     array length = 1 to map_tag_max
     set entire array to -1 as initial values
   for hash option:
     map_nhash = length of hash table
     map_nbucket = # of hash buckets, prime larger than map_nhash * 2
       so buckets will only be filled with 0 or 1 atoms on average
------------------------------------------------------------------------- */

void AtomKokkos::map_init(int check)
{
  // check for new map style if max atomID changed (check = 1 = default)
  // recreate = 1 if must delete old map and create new map
  // recreate = 0 if can re-use old map w/out realloc and just adjust settings
  // map_maxarray/map_nhash initially -1, to force recreate even when no atoms

  int recreate = 0;
  if (check) recreate = map_style_set();

  if (map_style == MAP_ARRAY && map_tag_max > map_maxarray)
    recreate = 1;
  else if (map_style == MAP_HASH && nlocal + nghost > map_nhash)
    recreate = 1;

  // if not recreating:
  // for array, initialize current map_tag_max values
  // for hash, set all buckets to empty, put all entries in free list

  if (!recreate) {
    if (lmp->kokkos->atom_map_classic) {
      if (map_style == MAP_ARRAY) {
        map_clear();
      } else {
        for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;
        map_nused = 0;
        map_free = 0;
        for (int i = 0; i < map_nhash; i++) map_hash[i].next = i + 1;
        if (map_nhash > 0) map_hash[map_nhash - 1].next = -1;
      }
    } else { 
      map_clear();
    }

    // recreating: delete old map and create new one for array or hash

  } else {
    map_delete();

    if (map_style == MAP_ARRAY) {
      map_maxarray = map_tag_max;
      memoryKK->create_kokkos(k_map_array, map_array, map_maxarray + 1, "atom:map_array");
      map_clear();
    } else {

      // map_nhash = max # of atoms that can be hashed on this proc
      // set to max of ave atoms/proc or atoms I can store
      // multiply by 2, require at least 1000
      // doubling means hash table will need to be re-init only rarely

      int nper = static_cast<int>(natoms / comm->nprocs);
      map_nhash = MAX(nper, nmax);
      map_nhash *= 2;
      map_nhash = MAX(map_nhash, 1000);

      if (lmp->kokkos->atom_map_classic) {
        // map_nbucket = prime just larger than map_nhash
        // next_prime() should be fast enough,
        //   about 10% of odd integers are prime above 1M

        map_nbucket = next_prime(map_nhash);

        // set all buckets to empty
        // set hash to map_nhash in length
        // put all hash entries in free list and point them to each other

        map_bucket = new int[map_nbucket];
        for (int i = 0; i < map_nbucket; i++) map_bucket[i] = -1;

        map_hash = new HashElem[map_nhash];
        map_nused = 0;
        map_free = 0;
        for (int i = 0; i < map_nhash; i++) map_hash[i].next = i + 1;
        map_hash[map_nhash - 1].next = -1;
      }

      k_map_hash = dual_hash_type(map_nhash);
    }
  }

  if (lmp->kokkos->atom_map_classic)
    if (map_style == MAP_ARRAY) k_map_array.modify_host();
}

/* ----------------------------------------------------------------------
   clear global -> local map for all of my own and ghost atoms
   for hash table option:
     global ID may not be in table if image atom was already cleared
------------------------------------------------------------------------- */

void AtomKokkos::map_clear()
{
  if (map_style == MAP_ARRAY) {
    if (lmp->kokkos->atom_map_classic) {
      Kokkos::deep_copy(k_map_array.h_view,-1);
      k_map_array.modify_host();
    } else {
      Kokkos::deep_copy(k_map_array.d_view,-1);
      k_map_array.modify_device();
    }
  } else {
    if (lmp->kokkos->atom_map_classic) {
      k_map_hash.h_view.clear();
      k_map_hash.modify_host();
      Atom::map_clear(); 
    } else {
      k_map_hash.d_view.clear();
      k_map_hash.modify_device();
    }
  }
}

/* ----------------------------------------------------------------------
   set global -> local map for all of my own and ghost atoms
   loop in reverse order so that nearby images take precedence over far ones
     and owned atoms take precedence over images
   this enables valid lookups of bond topology atoms
   for hash table option:
     if hash table too small, re-init
     global ID may already be in table if image atom was set
------------------------------------------------------------------------- */

void AtomKokkos::map_set()
{
  if (lmp->kokkos->atom_map_classic)
    map_set_host();
  else
    map_set_device();
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::map_set_device()
{
  int nall = nlocal + nghost;

  // possible reallocation of sametag must come before loop over atoms
  // since loop sets sametag

  if (nall > max_same) {
    max_same = nall + EXTRA;
    memoryKK->destroy_kokkos(k_sametag, sametag);
    memoryKK->create_kokkos(k_sametag, sametag, max_same, "atom:sametag");
  }

  if (map_style == MAP_HASH) {

    // if this proc has more atoms than hash table size, call map_init()
    //   call with 0 since max atomID in system has not changed
    // possible reallocation of sametag must come after map_init(),
    //   b/c map_init() may invoke map_delete(), whacking sametag

    if (nall > map_nhash) map_init(0);
  }

  atomKK->sync(Device, TAG_MASK);

  auto d_tag = atomKK->k_tag.d_view;
  auto d_sametag = k_sametag.d_view;

  // sort by tag

  int nmax = atom->nmax;

  int realloc_flag = 0;
  if (!d_tag_sorted.data() || (int)d_tag_sorted.extent(0) < nmax) {
    MemKK::realloc_kokkos(d_tag_sorted,"atom:tag_sorted",nmax);
    MemKK::realloc_kokkos(d_i_sorted,"atom:i_sorted",nmax);
    realloc_flag = 1;
  }

  h_tag_min() = MAXTAGINT;
  h_tag_max() = 0;

  Kokkos::deep_copy(d_tag_min_max,h_tag_min_max);

  auto l_tag_sorted = d_tag_sorted;
  auto l_i_sorted = d_i_sorted;
  auto l_tag_min = d_tag_min;
  auto l_tag_max = d_tag_max;
  int map_style_array = (map_style == MAP_ARRAY);

  Kokkos::parallel_for(nall, LAMMPS_LAMBDA(int i) {
    l_i_sorted(i) = i;
    tagint tag_i = d_tag(i);
    l_tag_sorted(i) = tag_i;
    Kokkos::atomic_min(&l_tag_min(),tag_i);
    Kokkos::atomic_max(&l_tag_max(),tag_i);
  });

  Kokkos::deep_copy(h_tag_min_max,d_tag_min_max);

  tagint min = h_tag_min();
  tagint max = h_tag_max();

  using MapKeyViewType = decltype(d_tag_sorted);
  using BinOpMap = Kokkos::BinOp1D<MapKeyViewType>;

  mapBinner = BinOpMap(nall, min, max);

  if (realloc_flag) {
    mapSorter = Kokkos::BinSort<MapKeyViewType, BinOpMap>(d_tag_sorted, 0, nall, mapBinner, true);
    MemKK::realloc_kokkos(mapSorter.bin_count_atomic,"Kokkos::SortImpl::BinSortFunctor::bin_count",nmax+1);
    Kokkos::deep_copy(mapSorter.bin_count_atomic,0);
    mapSorter.bin_count_const = mapSorter.bin_count_atomic;
    MemKK::realloc_kokkos(mapSorter.bin_offsets,"Kokkos::SortImpl::BinSortFunctor::bin_offsets",nmax+1);
    MemKK::realloc_kokkos(mapSorter.sort_order,"Kokkos::SortImpl::BinSortFunctor::sort_order",nmax);
  } else {
    Kokkos::deep_copy(mapSorter.bin_count_atomic,0);
    mapSorter.bin_op = mapBinner;
    mapSorter.range_begin = 0;
    mapSorter.range_end = nall;
  }

  mapSorter.create_permute_vector(LMPDeviceType());
  mapSorter.sort(LMPDeviceType(), d_tag_sorted, 0, nall);
  mapSorter.sort(LMPDeviceType(), d_i_sorted, 0, nall);

  auto d_map_array = k_map_array.d_view;
  auto d_map_hash = k_map_hash.d_view;
  d_map_hash.clear();

  auto d_error_flag = k_error_flag.d_view;
  Kokkos::deep_copy(d_error_flag,0);

  // for each tag find:
  //  neighboring atoms with closest local id for sametag
  //  atom with smallest local id for atom map

  Kokkos::parallel_for(nall, LAMMPS_LAMBDA(int ii) {
    const int i = l_i_sorted(ii);
    const tagint tag_i = l_tag_sorted(ii);

    int i_min = i;
    int i_closest = MAXSMALLINT;

    // search atoms with same tag in the forward direction

    int jj = ii+1;
    int closest_flag = 0;

    while (jj < nall) {
      const tagint tag_j = l_tag_sorted(jj);
      if (tag_j != tag_i) break;
      const int j = l_i_sorted(jj);
      i_min = MIN(i_min,j);
      if (j > i) {
        i_closest = MIN(i_closest,j);
        closest_flag = 1;
      }
      jj++;
    }

    // search atoms with same tag in the reverse direction

    jj = ii-1;

    while (jj >= 0) {
      const tagint tag_j = l_tag_sorted(jj);
      if (tag_j != tag_i) break;
      const int j = l_i_sorted(jj);
      i_min = MIN(i_min,j);
      if (j > i) {
        i_closest = MIN(i_closest,j);
        closest_flag = 1;
      }
      jj--;
    }

    if (!closest_flag)
      i_closest = -1;

    d_sametag(i) = i_closest;

    if (i == i_min) {
      if (map_style_array)
        d_map_array(tag_i) = i_min;
      else {
        auto insert_result = d_map_hash.insert(tag_i, i_min);
        if (insert_result.failed()) d_error_flag() = 1;
      }
    }

  });

  auto h_error_flag = k_error_flag.h_view;
  Kokkos::deep_copy(h_error_flag,d_error_flag);

  if (h_error_flag())
    error->one(FLERR,"Failed to insert into Kokkos hash atom map");

  k_sametag.modify_device();
  k_sametag.sync_host();

  if (map_style == MAP_ARRAY)
    k_map_array.modify_device();
  else
    k_map_hash.modify_device();
}

/* ---------------------------------------------------------------------- */

void AtomKokkos::map_set_host()
{
  int nall = nlocal + nghost;

  atomKK->sync(Host, TAG_MASK);
  k_sametag.sync_host();

  if (map_style == MAP_ARRAY) {
    k_map_array.sync_host();

    // possible reallocation of sametag must come before loop over atoms
    // since loop sets sametag

    if (nall > max_same) {
      max_same = nall + EXTRA;
      memoryKK->destroy_kokkos(k_sametag, sametag);
      memoryKK->create_kokkos(k_sametag, sametag, max_same, "atom:sametag");
    }

    for (int i = nall - 1; i >= 0; i--) {
      sametag[i] = map_array[tag[i]];
      map_array[tag[i]] = i;
    }

  } else {

    // if this proc has more atoms than hash table size, call map_init()
    //   call with 0 since max atomID in system has not changed
    // possible reallocation of sametag must come after map_init(),
    //   b/c map_init() may invoke map_delete(), whacking sametag

    if (nall > map_nhash) map_init(0);
    if (nall > max_same) {
      max_same = nall + EXTRA;
      memoryKK->destroy_kokkos(k_sametag, sametag);
      memoryKK->create_kokkos(k_sametag, sametag, max_same, "atom:sametag");
    }

    int previous, ibucket, index;
    tagint global;

    for (int i = nall - 1; i >= 0; i--) {
      sametag[i] = Atom::map_find_hash(tag[i]);

      // search for key
      // if found it, just overwrite local value with index

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index > -1) {
        map_hash[index].local = i;
        continue;
      }

      // take one entry from free list
      // add the new global/local pair as entry at end of bucket list
      // special logic if this entry is 1st in bucket

      index = map_free;
      map_free = map_hash[map_free].next;
      if (previous == -1)
        map_bucket[ibucket] = index;
      else
        map_hash[previous].next = index;
      map_hash[index].global = global;
      map_hash[index].local = i;
      map_hash[index].next = -1;
      map_nused++;
    }

    // Copy to Kokkos hash

    // use "view" template method to avoid unnecessary deep_copy

    auto h_map_hash = k_map_hash.view<LMPHostType>();
    h_map_hash.clear();

    for (int i = nall - 1; i >= 0; i--) {

      // search for key
      // if don't find it, done

      previous = -1;
      global = tag[i];
      ibucket = global % map_nbucket;
      index = map_bucket[ibucket];
      while (index > -1) {
        if (map_hash[index].global == global) break;
        previous = index;
        index = map_hash[index].next;
      }
      if (index == -1) continue;

      int local = map_hash[index].local;

      auto insert_result = h_map_hash.insert(global, local);
      if (insert_result.failed()) error->one(FLERR, "Kokkos::UnorderedMap insertion failed");
    }
  }

  k_sametag.modify_host();
  if (map_style == MAP_ARRAY)
    k_map_array.modify_host();
  else if (map_style == MAP_HASH)
    k_map_hash.modify_host();
}

/* ----------------------------------------------------------------------
   set global to local map for one atom
   for hash table option:
     global ID may already be in table if atom was already set
   called by Special class
------------------------------------------------------------------------- */

void AtomKokkos::map_one(tagint global, int local)
{
  if (lmp->kokkos->atom_map_classic)
    return Atom::map_one(global,local);

  if (map_style == MAP_ARRAY) {
    k_map_array.sync_host();
    k_map_array.h_view[global] = local;
  } else {
    k_map_hash.sync_host();
    auto& h_map_hash = k_map_hash.h_view;

    auto insert_result = h_map_hash.insert(global, local);
    if (insert_result.existing())
      h_map_hash.value_at(h_map_hash.find(global)) = local;
    else if (insert_result.failed())
      error->one(FLERR,"Failed to insert into Kokkos hash atom map");
  }
}

/* ----------------------------------------------------------------------
   lookup global ID in hash table, return local index
   called by map() in atom.h
------------------------------------------------------------------------- */

int AtomKokkos::map_find_hash(tagint global)
{
  k_map_hash.sync_host();
  auto& h_map_hash = k_map_hash.h_view;

  int local = -1;
  auto index = h_map_hash.find(global);
  if (h_map_hash.valid_at(index))
    local = h_map_hash.value_at(index);
  return local;
}

/* ----------------------------------------------------------------------
   free the array or hash table for global to local mapping
------------------------------------------------------------------------- */

void AtomKokkos::map_delete()
{
  memoryKK->destroy_kokkos(k_sametag, sametag);
  sametag = nullptr;

  if (map_style == MAP_ARRAY) {
    memoryKK->destroy_kokkos(k_map_array, map_array);
    map_array = nullptr;
  } else
    k_map_hash = dual_hash_type();

  if (lmp->kokkos->atom_map_classic)
    Atom::map_delete();
}
