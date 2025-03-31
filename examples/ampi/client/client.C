#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ccs-client.h"

#define SHRINK 1
#define EXPAND 0
#define OP EXPAND

#define BUF 255

int main (int argc, char **argv)
{
    int OLDNPROCS, NEWNPROCS;

    if (argc < 5) {
        printf("Usage: %s <hostname> <port> <oldprocs> <newprocs> \n", argv[0]);
        return 1;
    }

    // Create a CcsServer and connect to the given hostname and port
    CcsServer server;
    char host[BUF], *bitmap;
    int i, port, cmdLen, mode;

    sprintf(host, "%s", argv[1]);
    sscanf(argv[2], "%d", &port);
    sscanf(argv[3], "%d", &OLDNPROCS);
    sscanf(argv[4], "%d", &NEWNPROCS);

    if( NEWNPROCS > OLDNPROCS)
        mode = EXPAND;
    else if(OLDNPROCS > NEWNPROCS)
        mode = SHRINK;
    else{
        printf("Error: Old and new PE number is the same!\n");
        return 0;
    }
    printf("Connecting to server %s %d\n", host, port);
    CcsConnect(&server, host, port, NULL);
    printf("Connected to server\n");

    // define some variables to test
    int numNodes, numPes, nodeFirst, nodeSize;
    // get their values
    numNodes = CcsNumNodes(&server);
    numPes = CcsNumPes(&server);
    nodeFirst = CcsNodeFirst(&server, 0);
    nodeSize = CcsNodeSize(&server, 0);
    
    printf("num nodes: %d\n", numNodes);
    printf("num PEs: %d\n", numPes);
    printf("First node: %d\n", nodeFirst);
    printf("node size: %d\n", nodeSize);
    
    char *sendWhat = "client message";

    cmdLen = OLDNPROCS * sizeof(char) + sizeof(int) + sizeof(char);
    bitmap = (char *) malloc(cmdLen);

    if (mode == EXPAND) {
        printf("Sending expand command.\n");
        for (i = 0; i < OLDNPROCS; i++) {
            bitmap[i] = 1;
        }
    }
    else {
        printf("Sending shrink command.\n");
        for (i = 0; i < OLDNPROCS; i++) {
            if (i < NEWNPROCS)
                bitmap[i] = 1;
            else
                bitmap[i] = 0;
        }
    }
    memcpy(&bitmap[OLDNPROCS], &NEWNPROCS, sizeof(int));
    bitmap[OLDNPROCS+sizeof(int)] = '\0';
    // CcsSendRequest(&server, "set_bitmap", 0, cmdLen, bitmap);
    CcsSendRequest(&server, "check_shr_exp_", 0, strlen(sendWhat)+1, sendWhat);

    printf("Waiting for reply...\n" );
    // CcsRecvResponse(&server, cmdLen, bitmap , 180);
    CcsRecvResponse(&server, strlen(sendWhat)+1, sendWhat , 180);
    printf("Reply received.\n");

    return 0;
}
