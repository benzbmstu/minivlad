#!/usr/local/bin/python3.7
import sys
import socket
import logging
import logging.config


def main():
    log_setting = {
        "version": 1,
        "handlers": {
            "fileHandler": {
                "class":"logging.FileHandler",
                "formatter":"newFormatter",
                "filename":"pyserver.log"
            },
            "consoleHandler": {
                "class": "logging.StreamHandler",
                "level": "DEBUG",
                "formatter": "newFormatter"
            }
        },
        "loggers": {
            "ipm": {
                "handlers":["fileHandler", "consoleHandler"],
                "level":"DEBUG"
            }
        },
        "formatters": {
            "generalFormatter": {
                "format":"%(filename)s[LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s"
            },
            "testFormatter": {
                "format":"%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            },
            "newFormatter": {
                "format":"[%(name)s: %(levelname)-4s] %(asctime)s %(message)s [LINE:%(lineno)d]"
            }
        }
    }
    logging.config.dictConfig(log_setting)
    log = logging.getLogger('pyserver')
    log.info('preparing')

    listen_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    listen_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listen_socket.bind(('127.0.0.1', 8000))
    listen_socket.listen(1)

    while True:
        log.info('listening...')
        sys.stdout.write('listening...\n')
        client_connection, client_address = listen_socket.accept()
        sys.stdout.write('cl_addr: {}'.format(client_address))
        while True:
            request_link = client_connection.recv(1024)
            if not request_link:
                break
            log.info('get some link')
            sys.stdout.write('got some link\n')
            result = process_request(request_link.decode('utf-8').strip())
            print('result : {}'.format(result))
            resp = str(result).encode('utf-8') + b'\n'
            log.info(resp)
            client_connection.send(resp)



if __name__=='__main__':
    main()