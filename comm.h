#ifndef __COMM_H__
#define __COMM_H__

#define PORT	10000
#define BLOCKSIZE 1

int  saferead(int sd, unsigned char buf[], int len);

int  safewrite(int sd, unsigned char buf[], int len);

int createClient(char *ip, int port);

int comm(int sd, unsigned char dat);



#endif // __COMM_H__
