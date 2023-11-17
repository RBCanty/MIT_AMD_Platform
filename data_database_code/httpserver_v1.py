
from http.server import BaseHTTPRequestHandler, HTTPServer
from os import curdir, sep
from urllib import parse
import json
import pymongo
from add_data_database import add_data_mongo
from query_data_database import query_data_mongo
from add_reaction_database import add_reaction_mongo
from update_database_document import add_model_information, update_data_mongo, complete_data_campaign
import datetime
from pprint import pprint
import os
import logging
import re
import ssl
import traceback
import glob


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
        # Anytime a client issues a Post request to the HTTP Server we open a log file
        # If the log file does not exist we make a directory and then start
        # TODO: Ensure that logging is formatting correctly consistently
        current_time = datetime.datetime.now()
        log_path = log_path_base + current_time.strftime("%Y%m%d")
        if not os.path.exists(log_path):
            os.makedirs(log_path)
            print('Made a new directory: ' + log_path)
        LOG_FILENAME = log_path + '\\execution.log'
        logging.basicConfig(filename=LOG_FILENAME, level=logging.INFO)
        
        # Checks to ensure that the request comes from a valid client ip address
        # Basically the 18.0.0.0 - 18.31.255.255 are MIT addresses for most of campus
        # TODO: Check to make sure that the CSAIL folks have access as well
        #if server_testing != True:
        #    valid_ip_addresses = '18\.[0-3]\d\.\d{1,3}\.\d{1,3}'
        #    if re.match(valid_ip_addresses, self.client_address[0]) == None:
        #        logging.info("POST request IP failed from,\nClient: %s", str(self.client_address))
        #        self.send_response(401)
        #        self.end_headers
        #        return
        
        # Confirms that the header from the request matches some "content-type"
        # Can be set to an arbitrary mutual client-server value
        # If the content type matches it converts the request to a dictionary
        # If not we return a 401 (Unauthorized) response to the client
        if self.headers.get('content-type') == 'application/authjson':
            length = int(self.headers.get('content-length'))
            message = json.loads(self.rfile.read(length))
            logging.info("POST request,\nClient: %s\nHeaders:\n%s \n Body:\n%s",\
                         str(self.client_address), str(self.headers).replace('\n',' <> '), \
                         str(json.dumps(message['request'])))
            
            # Try to handle the message contents
            # Starts by attempting to authenticate user on the Mongo Database
            # Returns a 403 (Forbidden) response on authentication failure
            try:
                myclient = pymongo.MongoClient('mongodb://' + message['auth']['user'] + ':' + message['auth']['password'] + \
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
            # TODO: Implement additional request types
            try:
                current_time = datetime.datetime.now()
                requesttype = message['request']['request_type']
                if requesttype == 'add_data':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    data_response = add_data_mongo(message['request'],myclient)
                    response_code = 200
                elif requesttype == 'query_data':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    data_response = query_data_mongo(message['request'],myclient)
                    response_code = 200
                elif requesttype == 'add_reaction':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    data_response = add_reaction_mongo(message['request'],myclient)
                    response_code = 200
                elif requesttype == 'add_model_information':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    data_response = add_model_information(message['request'],myclient)
                    response_code = 200
                elif requesttype == 'update_model_information':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    data_response = update_data_mongo(message['request'],myclient)
                    response_code = 200
                elif requesttype == 'complete_data_campaign':
                    print('Request type known: ' + requesttype + ' ' + current_time.strftime("%m%d%Y_%H%M%S"))
                    data_response = complete_data_campaign(message['request'],myclient)
                    response_code = 200
                # TODO: Add key-values for model_information
                # What is needed for things to work?
                # Needs to query the collection to look for new jobs
                # Needs to send jobs to the client
                # Needs to allow job documents to be updated
                else:
                    data_response = ['Error', {'print_statements': ['Request not recognized: ' + message['request']['request_type']]}]
                    response_code = 405
            except KeyError:
                print(traceback.format_exc())
                data_response = ['Error', {'print_statements': ['Data request not properly formatted']}]
                response_code = 400
            backup_files = glob.glob('.\backups\*.json')
            for backup_file in backup_files:
                if 'database_backup' not in backup_file:
                    continue
                time_elapsed = datetime.strptime(backup_file.split('_')[-1].split('.')[0], '%Y%m%d') - current_time
                if time_elapsed.total_seconds < 24 * 60 * 60:
                    break
            else:
                backup_filename = r'.\backups\database_backup_%s.json' % current_time.strftime("%Y%m%d")
                json_backup = {}
                collections = myclient.amd_database.command('listCollections')
                for collection in collections['cursor']['firstBatch']:
                    json_backup[collection['name']] = {}
                    document_cursor = myclient.amd_database[collection['name']].find({})
                    for document in document_cursor:
                        json_backup[collection['name']][str(document['_id'])] = document
                        json_backup[collection['name']][str(document['_id'])]['_id'] = str(document['_id'])
                
                with open(backup_filename, 'w') as jsonfile:
                    json.dump(json_backup, jsonfile)
            myclient.close()
        else:
            logging.info("Broken POST request,\nPath: %s\nHeaders:\n%s\n", str(self.path), str(self.headers))
            self.send_response(401)
            self.end_headers()
            return
        
        # If the code makes it to this point, something has happened on the database
        # We will take the database response and add the original request to it
        # Then we get everything together and send back a response to the client
        data_response[1]['server_request'] = message['request']
        logging.info(str(response_code) + '\n')
        logging.shutdown()
        self.send_response(response_code)
        self.send_header('Content-type', 'bfd1f418cbc85ac1d4f9fff50e291003')
        self.end_headers()
        # print(data_response)
        self.wfile.write(json.dumps(data_response).encode(encoding='utf_8'))
        # At this point the request handler is done
 
try:
    # Set the location of the log files, if base path does not exist all requests
    # will terminate with error responses and clients get forcibly disconnect
    # responses from the server
    log_path_base = 'C:\\Database\\logs\\'
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
    # Deafult FTIR Values : 01.23.45.678:20123
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
    print('Started httpserver on: %s:%s' % (HTTP_ADDRESS, str(PORT_NUMBER)))
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


