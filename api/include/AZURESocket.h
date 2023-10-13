#ifndef AZURESOCKET_H
#define AZURESOCKET_H

#include <iostream>
#include <string>
#include <cstring>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>

#include "AZUREAPI.h"

const int BUFFER_SIZE = 10000;

class AZURESocket {

private:
    int port_;
    int serverSocket_;
    int clientSocket_;
    struct sockaddr_in serverAddress_;
    struct sockaddr_in clientAddress_;

    AZUREAPI* api_;

public:
    AZURESocket(int port, AZUREAPI* api): port_(port), api_( api ) {
      serverSocket_ = -1;
      clientSocket_ = -1;
    };

    ~AZURESocket(){
      close(serverSocket_);
      close(clientSocket_);
    };

    bool start( );

    bool sendPacket( vector_r response );
    bool sendPacket( std::string );
    bool sendPacket( std::vector<bool> response );

};


#endif