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
    #logging.config.dictConfig(log_setting)
    #log = logging.getLogger('pyserver')
    #log.info('preparing')

    listen_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    listen_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listen_socket.bind(('127.0.0.1', 8081))
    listen_socket.listen(1)

    while True:
        #log.info('listening...')
        sys.stdout.write('listening...\n')
        client_connection, client_address = listen_socket.accept()

        # проверка установки соединения (если не удалось, то continue)
        sys.stdout.write('cl_addr: {}\n'.format(client_address))

        #client_connection.send('Hello world!'.encode('utf-8'))
        #client_connection.close() # close connection
        #continue

        while True:
            request_link = client_connection.recv(1024) # ожидаю данные от клиента именно 1024 байт
            #if not request_link:
            #    break
            # log.info('get some link')
            sys.stdout.write('got some link\n')
            result = request_link.decode('utf-8').strip()
            print('result : {}'.format(result))
            #resp = str(result).encode('utf-8') + b'\n'
            #print("resp", resp)
            #log.info(resp)
            #client_connection.send(resp)
            client_connection.send('Hello world!'.encode('utf-8'))
            client_connection.close()
            break


def part_2():
    my_socket = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
    my_socket.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
    my_socket.bind(('127.0.0.1', 8081))
    my_socket.listen(1)

    #socket creation and other functions done here (see previous post)
    while True:
        sys.stdout.write('listening...\n')
        connection,address = my_socket.accept()
        request = connection.recv(1024).decode('utf-8')
        string_list = request.split(' ')     # Split request from spaces

        method = string_list[0] # First string is a method
        requesting_file = string_list[1] #Second string is request file
        print('Client request ',requesting_file)


        myfile = requesting_file.split('?')[0] # After the "?" symbol not relevent here       
        myfile = myfile.lstrip('/')
        if(myfile == '/'):
            myfile = 'index.html'    # Load index file as default


        try:
            file = open(myfile,'rb') # open file , r => read , b => byte format
            response = file.read()
            file.close()
 
            header = 'HTTP/1.1 200 OK\n'
 
            if(myfile.endswith(".jpg")):
                mimetype = 'image/jpg'
            elif(myfile.endswith(".css")):
                mimetype = 'text/css'
            else:
                mimetype = 'text/html'
 
            header += 'Content-Type: '+str(mimetype)+'<strong>\n\n</strong>'
 
        except Exception as e:
            header = 'HTTP/1.1 404 Not Found\n\n'
            response = '<html>\
                          <body>\
                            <center>\
                             <h3>Error 404: File not found</h3>\
                             <p>Python HTTP Server</p>\
                            </center>\
                          </body>\
                        </html>'.encode('utf-8')

        final_response = header.encode('utf-8')
        final_response += response
        connection.send(final_response)
        connection.close()


def part_3():
    listen_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    listen_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listen_socket.bind(('127.0.0.1', 8081))
    listen_socket.listen(1)

    while True:
        sys.stdout.write('listening...\n')
        client_connection, client_address = listen_socket.accept()

        # проверка установки соединения (если не удалось, то continue)
        sys.stdout.write('cl_addr: {}\n'.format(client_address))

        while True:
            request_link = client_connection.recv(1024) # ожидаю данные от клиента именно 1024 байт
            if not request_link:
                break
            sys.stdout.write('got some link\n')
            result = request_link.decode('utf-8').strip()
            print('result : {}'.format(result))

            header = 'HTTP/1.1 200 OK\n\n'
            response = '<html>\
                          <body>\
                            <center>\
                             <p>Python HTTP Server</p>\
                            </center>\
                          </body>\
                        </html>'.encode('utf-8')

            final_response = header.encode('utf-8')
            final_response += response
            client_connection.send(final_response)
            client_connection.close()

            break



def set_proc_name(newname):
    from ctypes import cdll, byref, create_string_buffer
    libc = cdll.LoadLibrary('libc.so.6')
    buff = create_string_buffer(len(newname)+1)
    buff.value = newname
    libc.prctl(15, byref(buff), 0, 0, 0)

def get_proc_name():
    from ctypes import cdll, byref, create_string_buffer
    libc = cdll.LoadLibrary('libc.so.6')
    buff = create_string_buffer(128)
    # 16 == PR_GET_NAME from <linux/prctl.h>
    libc.prctl(16, byref(buff), 0, 0, 0)
    return buff.value

if __name__=='__main__':
    # sys.argv[0] == 'python'
    # outputs 'python'
    #get_proc_name()
    #set_proc_name('testing yeah')
    #sys.argv[0] == 'testing_yeah'
    # outputs 'testing yeah'

    #get_proc_name()
    part_3()



