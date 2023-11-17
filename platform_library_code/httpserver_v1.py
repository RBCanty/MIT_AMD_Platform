# -*- coding: utf-8 -*-
"""
Created on Mon, July 26th 2021

HTTP Server for interfacing with the Platform Mongo Database

@author: Brent
"""
from http.server import BaseHTTPRequestHandler, HTTPServer
from os import curdir, sep
from urllib import parse
import json, pymongo, datetime, os, logging, re, yaml
from platform_library_functions import without_keys, update_database_field, add_to_database, query_database, query_collection, query_document, move_database_document, update_database_document, query_platform_location
from pprint import pprint
from copy import deepcopy

# This class will handle incoming requests from the client
# TODO: Add a SSL handshake authentication to the socket
# TODO: Would starting a persistent connection be better? 
class myHandler(BaseHTTPRequestHandler):
    def _set_headers(self):
        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        
    def do_GET(self):
        # TODO: Implement get requests to return data? Possibly?
        # Right now the get requests are just a stub for future coding
        if server_testing != True:
            current_time = datetime.datetime.now()
            log_path = log_path_base + current_time.strftime("%Y%m%d")
            if not os.path.exists(log_path):
                os.makedirs(log_path)
                print('Made a new directory: ' + log_path)
            LOG_FILENAME = log_path + '\\execution.log'
            logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
            logging.info("GET request,\nPath: %s\nHeaders:\n%s\n", str(self.path), str(self.headers))
            logging.shutdown()
            self.send_response(400)
            self.end_headers()
            self.wfile.write(b'Hello this is just a stub right now!')
        else:
            self.send_response(400)
            self.end_headers
            return

    def do_POST(self):
        if self.headers.get('content-type') == 'application/authjson':
            length = int(self.headers.get('content-length'))
            message = json.loads(self.rfile.read(length))
            logging.info("POST request,\nClient: %s\nHeaders:\n%s \n Body:\n%s",\
                         str(self.client_address), str(self.headers).replace('\n',' <> '), \
                         str(json.dumps(message['request'])))
            try:
                platform_mongodb = pymongo.MongoClient('mongodb://' + message['auth']['user'] + ':' + message['auth']['password'] + \
                                           '@localhost:' + message['auth']['port'] + '/admin' + '?maxIdleTimeMS=100000')
                logging.info(message['auth']['user'] + ' connected to MongoDB!')
            except:
                logging.info(message['auth']['user'] + 'Failed connection to the MongoDB...')
                self.send_response(403)
                self.end_headers()
                return
            # Check to make sure that there is at least a request type in the request
            # With no request type present we return a 400 (Bad Request) response
            # If there is we will try to handle the request type
            # If not we return a 405 (Method Not Allowed) response
            data_response = {}
            try:
                current_time = datetime.datetime.now()
                requesttype = message['request']['request_type']
                if requesttype == 'add_to_database':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    config_filepath = r'C:\Users\admin\Documents\AMD_Control_Platform-main\venv\Shell_AH\Evoware_API\config.yaml'
                    with open(config_filepath, 'r') as yamlfile:
                        config = yaml.load(yamlfile, Loader = yaml.Loader)
                    carriers_labware_filepath = config['paths/files']['carriers_labware_path']
                    with open(carriers_labware_filepath, 'r') as jsonfile:
                        carriers_labware = json.load(jsonfile)
                    internal_request = deepcopy(message['request']['incoming_dict'])
                    request_response = add_to_database(internal_request, message['request']['collection'],
                                                       carriers_labware, platform_mongodb)
                    response_code = 200
                elif requesttype == 'query_collection':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    request_response = query_collection(message['request']['collection'], platform_mongodb)
                    response_code = 200
                elif requesttype == 'query_document':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    request_response = query_document(message['request']['collection'], message['request']['search_term'],
                                                      message['request']['search_field'], platform_mongodb)
                    response_code = 200
                elif requesttype == 'move_database_document':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    request_response = move_database_document(message['request']['doc_id'], message['request']['current_collection'],
                                                              message['request']['new_collection'], platform_mongodb)
                    response_code = 200
                elif requesttype == 'update_database_document':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    request_response = update_database_document(message['request']['previous_document'], message['request']['new_document'],
                                                              message['request']['collection'], platform_mongodb)
                    pprint(request_response)
                    response_code = 200
                    logging.info(str(request_response))
                elif requesttype == 'query_platform_location':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    request_response = query_platform_location(message['request']['location_tag'], message['request']['sublocation'],
                                                               platform_mongodb)
                    response_code = 200
                elif requesttype == 'update_database_field':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    request_response = update_database_field(message['request']['document_name'], message['request']['update_field'],
                                                             message['request']['old_value'], message['request']['new_value'],
                                                             platform_mongodb)
                    response_code = 200
                else:
                    data_response['data_errors'] = 'Request not recognized: ' + message['request']['request_type']
                    request_response = ['Request not recogized: %s' % message['request']['request_type']]
                    response_code = 405
            except Exception as e:
                error_path = error_path_base + current_time.strftime("%Y%m%d")
                if not os.path.exists(error_path):
                    os.makedirs(error_path)
                    print('Made a new directory: ' + error_path)
                error_filepath = error_path + '\\error_log.txt'
                with open(error_filepath, 'w') as error_file:
                    error_file.write(r'<>----- ¯\_(ツ)_/¯ -----<>\n')
                    error_file.write(str(e))
                    error_file.write('\n')
                    error_file.write(r'<>----- ¯\_(ツ)_/¯ -----<>\n')
                data_response = {}
                data_response['data_errors'] = 'Data request not properly formatted'
                request_response = ['Data request not properly formatted']
                response_code = 400
            platform_mongodb.close()
        else:
            logging.info("Broken POST request,\nPath: %s\nHeaders:\n%s\n", str(self.path), str(self.headers))
            self.send_response(401)
            self.end_headers()
            return
        data_response['server_request'] = message['request']
        data_response['request_return'] = request_response
        logging.info(str(response_code) + '\n')
        logging.shutdown()
        self.send_response(response_code)
        self.send_header('Content-type', 'application/json')
        self.end_headers()
        self.wfile.write(json.dumps(data_response).encode(encoding='utf_8'))
        # At this point the request handler is done
 
try:
    # Set the location of the log files, if base path does not exist all requests
    # will terminate with error responses and clients get forcibly disconnect
    # responses from the server
    log_path_base = r'C:\Users\admin\Documents\amd_project_supporting_docs\platform_library\http_logs\\'
    error_path_base = r'C:\Users\admin\Documents\amd_project_supporting_docs\platform_library\error_logs\\'
    current_time = datetime.datetime.now()
    log_path = log_path_base + current_time.strftime("%Y%m%d")
    if not os.path.exists(log_path):
        os.makedirs(log_path)
        print('Made a new directory: ' + log_path)
    LOG_FILENAME = log_path + '\\execution.log'
    logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
    
    # Just to switch between testing on local machines to network accessible
    # this changes the Port Number and HTTP Address to connect to, this also
    # changes the IP address requirements in the POST handler and opens the 
    # GET handler process tree as well
    # Default Testing Values : 127.0.0.1:9999
    # Deafult FTIR Values : 12.34.56.789:20123
    server_testing = False
    if server_testing == True:
            PORT_NUMBER = 9999
            HTTP_ADDRESS = '127.0.0.1'
    else:
        PORT_NUMBER = 20123
        HTTP_ADDRESS = '127.0.0.1'
    
    # Now we start the HTTP Server and have it serve forever
    # We just open and close the log files as we need them to declutter memory
    # TODO: Implement an SSL handshake on client and server for encryption
    httpd = HTTPServer((HTTP_ADDRESS, PORT_NUMBER), myHandler)
    #httpd.socket = ssl.wrap_socket(httpd.socket, keyfile = 'E:\\Version 1 AMD Database\\key1.pem', certfile='E:\\Version 1 AMD Database\\cert1.pem')
    print('Started httpserver on %s:%s' % (str(HTTP_ADDRESS), str(PORT_NUMBER)))
    logging.info('Starting httpd...\n')
    logging.shutdown()
    httpd.serve_forever()    

except KeyboardInterrupt:
    # When the server recieves a ^C command, it will stop serving and shutdown
    # otherwise as long as the terminal is open the server will continue even 
    # if there is an error raised
    print('^C received, shutting down the web server')
    httpd.socket.close()
    current_time = datetime.datetime.now()
    log_path = log_path_base + current_time.strftime("%Y%m%d")
    if not os.path.exists(log_path):
        os.makedirs(log_path)
        print('Made a new directory: ' + log_path)
    LOG_FILENAME = log_path + '\\execution.log'
    logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
    logging.info('Stopping httpd...\n')
