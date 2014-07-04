/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "fix_ipi.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "kspace.h"
#include "modify.h"
#include "compute.h"
#include "neighbor.h"
#include "domain.h"
#include "compute_pressure.h"
#include "errno.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// socket interface
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <netdb.h>


#define MSGLEN 12
namespace sockets {

void error(const char *msg)
// Prints an error message and then exits.

{   perror(msg);  exit(0);   }

void open_socket(int *psockfd, int* inet, int* port, char* host)
/* Opens a socket.

Note that fortran passes an extra argument for the string length, but this is
ignored here for C compatibility.

Args:
   psockfd: The id of the socket that will be created.
   inet: An integer that determines whether the socket will be an inet or unix
      domain socket. Gives unix if 0, inet otherwise.
   port: The port number for the socket to be created. Low numbers are often
      reserved for important channels, so use of numbers of 4 or more digits is
      recommended.
   host: The name of the host server.
*/

{
   int sockfd, portno, n;
   struct hostent *server;

   fprintf(stderr, "Connection requested %s, %d, %d\n", host, *port, *inet);
   struct sockaddr * psock; int ssock;
   struct sockaddr_in serv_addr_in; struct sockaddr_un serv_addr_un;
   if (*inet>0)
   {
      psock=(struct sockaddr *)&serv_addr_in;     ssock=sizeof(serv_addr_in);
      sockfd = socket(AF_INET, SOCK_STREAM, 0);
      if (sockfd < 0)  error("ERROR opening socket");

      server = gethostbyname(host);
      if (server == NULL)
      {
         fprintf(stderr, "ERROR, no such host %s \n", host);
         exit(0);
      }

      bzero((char *) &serv_addr_in, ssock);
      serv_addr_in.sin_family = AF_INET;
      bcopy((char *)server->h_addr, (char *)&serv_addr_in.sin_addr.s_addr, server->h_length);
      serv_addr_in.sin_port = htons(*port);
   }
   else
   {
      psock=(struct sockaddr *)&serv_addr_un;     ssock=sizeof(serv_addr_un);
      sockfd = socket(AF_UNIX, SOCK_STREAM, 0);
      bzero((char *) &serv_addr_un, ssock);
      serv_addr_un.sun_family = AF_UNIX;
      strcpy(serv_addr_un.sun_path, "/tmp/ipi_");
      strcpy(serv_addr_un.sun_path+9, host);
   }
   if (connect(sockfd, psock, ssock) < 0) error("ERROR connecting");

   *psockfd=sockfd;
}

void writebuffer(int sockfd, char *data, int len)
/* Writes to a socket.

Args:
   psockfd: The id of the socket that will be written to.
   data: The data to be written to the socket.
   plen: The length of the data in bytes.
*/
{
   int n;

   n = write(sockfd,data,len);
   if (n < 0) error("ERROR writing to socket");
}


void readbuffer(int sockfd, char *data, int len)
/* Reads from a socket.

Args:
   psockfd: The id of the socket that will be read from.
   data: The storage array for data read from the socket.
   plen: The length of the data in bytes.
*/

{
   int n, nr;


   n = nr = read(sockfd,data,len);

   while (nr>0 && n<len )
   {  nr=read(sockfd,&data[n],len-n); n+=nr; }

   if (n == 0) error("ERROR reading from socket");
}


}; // ends namespace socket

/* ---------------------------------------------------------------------- */

FixIPI::FixIPI(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  /* format for fix:
   *  fix  num  group_id  driver  host port [unix]
   */
  if (strcmp(style,"ipi") != 0 && narg < 5)
    error->all(FLERR,"Illegal fix ipi command");

  //box_change = 1;
  //box_change_size = 1;

  strcpy(host, arg[3]);
  port=atoi(arg[4]);

  if (narg > 5 && strcmp(arg[5],"unix") ==0 ) inet=0; else inet=1;

  MPI_Comm_rank(world,&me);
  if (me==0) master=1; else master=0;

  hasdata=0; bsize=0;


  char** newarg = new char*[3];
  newarg[0] = (char *) "ipi_temp";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
//  tflag = 1;

  newarg = new char*[5];
  newarg[0] = (char *) "ipi_press";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = (char *) "ipi_temp";
  newarg[4] = (char *) "virial";
  modify->add_compute(5,newarg);
  delete [] newarg;
//  pflag = 1;
}

/* ---------------------------------------------------------------------- */

int FixIPI::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIPI::init()
{
  if (master) sockets::open_socket(&ipisock, &inet, &port, host);
  else ipisock=0;
  //! should check for success in socket opening

  // asks for evaluation of PE at first step
  modify->compute[modify->find_compute("thermo_pe")]->invoked_scalar = -1;
  modify->addstep_compute_all(update->ntimestep + 1);

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // makes sure that neighbor lists are re-built at each step (cannot make assumptions when cycling over beads!)
  neighbor->delay = 0;
  neighbor->every = 1;
}

void FixIPI::initial_integrate(int vflag)
{
  char header[MSGLEN+1];

  if (hasdata)
    error->all(FLERR, "Driver got out of sync in initial_integrate and will die!");

  double cellh[9], cellih[9];
  int nat;
  if (master)   {

    while (true) {
      sockets::readbuffer(ipisock, header, MSGLEN); header[MSGLEN]=0;

      if (strcmp(header,"STATUS      ") == 0 )
        sockets::writebuffer(ipisock,"READY       ",MSGLEN);
      else break;
    }

    if (strcmp(header,"POSDATA     ") == 0 )  {
      sockets::readbuffer(ipisock, (char*) cellh, 9*sizeof(double));
      sockets::readbuffer(ipisock, (char*) cellih, 9*sizeof(double));
      sockets::readbuffer(ipisock, (char*) &nat, sizeof(int));

      if (bsize==0) { bsize=3*nat; buffer = new double[bsize]; }
      else if (bsize!=3*nat)
        error->all(FLERR, "Number of atoms changed along the way.");

      sockets::readbuffer(ipisock, (char*) buffer, sizeof(double)*bsize);
    }
    else error->all(FLERR, "Wrapper did not send positions, I will now die!");

  }

  //shares the atomic coordinates with everyone
  //MPI_Bcast(&nat,1,MPI_INTEGER,0,world);
  MPI_Bcast(&nat,1,MPI_INT,0,world);
  //must also allocate the buffer on the non-head nodes
  if (bsize==0) { bsize=3*nat; buffer = new double[bsize]; }
  MPI_Bcast(cellh,9,MPI_DOUBLE,0,world);
  MPI_Bcast(cellih,9,MPI_DOUBLE,0,world);
  MPI_Bcast(buffer,bsize,MPI_DOUBLE,0,world);


  //updates atomic coordinates based on the data received
  //!! SHOULD ALSO WORK OUT HOW TO UPDATE CELL!
  double *boxhi = domain->boxhi;
  double *boxlo = domain->boxlo;
  double posconv;
  posconv=0.52917721*force->angstrom;
  boxlo[0] = 0;
  boxlo[1] = 0;
  boxlo[2] = 0;
  boxhi[0] = cellh[0]*posconv;
  boxhi[1] = cellh[4]*posconv;
  boxhi[2] = cellh[8]*posconv;
  domain->xy = cellh[1]*posconv;
  domain->xz = cellh[2]*posconv;
  domain->yz = cellh[5]*posconv;
  domain->reset_box();
  domain->box_change = 1;
  if (kspace_flag) force->kspace->setup();
  //pressure->addstep(update->ntimestep+1);

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < nlocal; i++) {
     if (mask[i] & groupbit) {
       x[i][0]=buffer[3*(atom->tag[i]-1)+0]*posconv;
       x[i][1]=buffer[3*(atom->tag[i]-1)+1]*posconv;
       x[i][2]=buffer[3*(atom->tag[i]-1)+2]*posconv;
     }
  }

  // compute PE. makes sure that it will be evaluated at next step
  modify->compute[modify->find_compute("thermo_pe")]->invoked_scalar = -1;
  modify->addstep_compute_all(update->ntimestep+1);

  hasdata=1;
}

/* ---------------------------------------------------------------------- */

void FixIPI::final_integrate()
{
  char header[MSGLEN+1];
  double vir[9], pot=0.0;
  double forceconv, potconv, posconv, pressconv, posconv3;
  potconv=3.1668152e-06/force->boltz;
  posconv=0.52917721*force->angstrom;
  posconv3=posconv*posconv*posconv;
  forceconv=potconv*posconv;
  pressconv=1/force->nktv2p;

  // compute for potential energy
  pot=modify->compute[modify->find_compute("thermo_pe")]->compute_scalar();
  pot*=potconv;
  //!should also get virial (configurational bit only!)

  if (!hasdata)
    error->all(FLERR, "Driver got out of sync in final_integrate and will die!");

  int nat=bsize/3;
  double **f= atom->f, lbuf[bsize];

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  for (int i = 0; i < bsize; ++i) lbuf[i]=0.0;
  for (int i = 0; i < nlocal; i++) {        // do these also contain ghost atoms???
     lbuf[3*(atom->tag[i]-1)+0]=f[i][0]*forceconv;  //must check units
     lbuf[3*(atom->tag[i]-1)+1]=f[i][1]*forceconv;
     lbuf[3*(atom->tag[i]-1)+2]=f[i][2]*forceconv;
  }
  MPI_Allreduce(lbuf,buffer,bsize,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0; i < 9; ++i) vir[i]=0.0;
  //Need the following lines in the input:
  int press_id = modify->find_compute("ipi_press");
  Compute* comp_p = modify->compute[press_id];
  comp_p->compute_vector();
  double myvol = domain->xprd*domain->yprd*domain->zprd/posconv3;
  //std::cerr<<"vol "<<myvol<<"\n";
  //for (int i = 0; i < 6; i++)
  vir[0] = comp_p->vector[0]*pressconv*myvol;
  vir[4] = comp_p->vector[1]*pressconv*myvol;
  vir[8] = comp_p->vector[2]*pressconv*myvol;
  vir[5] = comp_p->vector[3]*pressconv*myvol;
  vir[2] = comp_p->vector[4]*pressconv*myvol;
  vir[1] = comp_p->vector[5]*pressconv*myvol;

  if (master) {
    while (true) {
      sockets::readbuffer(ipisock, header, MSGLEN); header[MSGLEN]=0;

      if (strcmp(header,"STATUS      ") == 0 )
        sockets::writebuffer(ipisock,"HAVEDATA    ",MSGLEN);
      else break;
    }

    if (strcmp(header,"GETFORCE    ") == 0 )  {

      sockets::writebuffer(ipisock,"FORCEREADY  ",MSGLEN);
      sockets::writebuffer(ipisock,(char*) &pot,sizeof(double));
      sockets::writebuffer(ipisock,(char*) &nat,sizeof(int));
      sockets::writebuffer(ipisock,(char*) buffer, bsize*sizeof(double));
      sockets::writebuffer(ipisock,(char*) vir,9*sizeof(double));
      nat=0;  sockets::writebuffer(ipisock,(char*) &nat,sizeof(int));
    }
    else
      error->all(FLERR, "Wrapper did not ask for forces, I will now die!");

  }

  hasdata=0;
}


