import requests
from pprint import pprint
import glob
import threading
import subprocess
import sys
import queue
from http.server import BaseHTTPRequestHandler, HTTPServer
import datetime
import json
import time
import os


def download_model_from_dropbox():
    # Implement a method to download files from separate HPC cluster
    return ['Error', 'Need to add code to handle this']


def prepare_local_retrain_dataset(model_type, model_version):
    model_directory = r'.\chemprop_models\%s' % model_type
    datafile_location = r'%s\data' % model_directory
    current_datafiles = glob.glob(r'%s\*.csv' % datafile_location)
    retraining_dataset = {}
    for current_datafile in current_datafiles:
        print('Working on file: %s' % current_datafile)
        with open(current_datafile, 'r') as infile:
            line_index = 0
            header_line = []
            for line in infile:
                if line_index == 0:
                    header_line = line
                    line_index += 1
                    continue
                if len(line.strip()) < 5:
                    continue
                elif 'smiles' in line.strip().lower():
                    continue
                line_split = line.strip().split(',')
                if line_split[0] not in retraining_dataset.keys():
                    retraining_dataset[line_split[0]] = []
                if line_split in retraining_dataset[line_split[0]]:
                    continue
                retraining_dataset[line_split[0]].append(line_split)

    retraining_dataset_filepath = r'%s\%s_dataset.csv' % (datafile_location, model_version)
    with open(retraining_dataset_filepath, 'w') as outfile:
        outfile.write(header_line)
        for item in retraining_dataset.keys():
            if len(retraining_dataset[item]) == 0:
                continue
            else:
                for data_point in retraining_dataset[item]:
                    outfile.write(','.join(data_point) + '\n')

    return ['Success', '', retraining_dataset_filepath]


def chemprop_thread_function(chemprop_command, job_name):
    print('\tReceived command: %s' % chemprop_command)
    process = subprocess.Popen(chemprop_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        print(line.decode(), end='')
    print('Finishing!: %s' % chemprop_command)


def retrain_model_from_new_data(model_type, model_version, save_location, new_datapoints=None):
    # Request new datapoints for model_type
    # Add the dataset to the current model dataset (check for duplicates)
    # Write the full dataset to a new filepath
    # Train chemprop model with the file name provided + _training
    # Change the model directory name to the model version name

    if new_datapoints is None:
        new_datapoints = {}
    if save_location == 'local':
        return_statement, return_object, dataset_filepath = prepare_local_retrain_dataset(model_type, model_version)
        if return_statement != 'Success':
            return ['Error', return_object]

        model_directory = os.path.abspath(r'.\chemprop_models\%s' % model_type)
        if model_type == 'photodeg_rf':
            chemprop_input_arguments = {'function': 'bash',
                                        'arguments': [r'%s\retrain_model_pipeline.sh' % model_directory,
                                                      r'%s\%s_training' % (model_directory, model_version),
                                                      os.path.abspath(dataset_filepath)]}

        else:
            return ['Error', 'Model retraining details for "%s" and "%s" not defined' % (save_location, model_type)]

        if chemprop_input_arguments['function'] == 'bash':
            chemprop_arguments = chemprop_input_arguments['arguments']
            chemprop_input_command = 'bash ' + ' '.join(chemprop_arguments)
        else:
            return ['Error', 'No implemented correctly for function type: %s' % chemprop_input_arguments['function']]
        chemprop_thread = threading.Thread(target=chemprop_thread_function,
                                           args=([chemprop_input_command, model_type]))
        chemprop_thread.start()
        chemprop_thread.join()

    return ['Success', '']


def retraining_thread(retraining_queue):
    time.sleep(10)
    while True:
        current_queue = retraining_queue.queue
        if len(current_queue[0]['queue']) == 0:
            time.sleep(60)
        else:
            current_queue = retraining_queue.get()
            current_item = current_queue['queue'][0]
            current_queue['running'].append(current_item)
            del(current_queue['queue'][0])
            retraining_queue.put(current_queue)
            print('\tStarting retraining on %s' % current_item['model_type'])
            match_items = ['model_type', 'model_version', 'retrain_location']
            print([current_item[key] for key in match_items])

            return_statement, return_object = retrain_model_from_new_data(current_item['model_type'],
                                                                          current_item['model_version'],
                                                                          current_item['retrain_location'])
            if return_statement != 'Success':
                sys.exit()
            current_queue = retraining_queue.get()
            current_queue['running'] = []
            retraining_queue.put(current_queue)
            time.sleep(60)


def MakeHandlerClass(retraining_queue):
    class RetrainingHandler(BaseHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            self.model_retraining_queue = retraining_queue
            super(RetrainingHandler, self).__init__(*args, **kwargs)

        def _set_headers(self):
            self.send_response(200)
            self.send_header('content-type', 'application/json')
            self.end_headers()

        def do_POST(self):
            length = int(self.headers.get('content-length'))
            message = json.loads(self.rfile.read(length))
            pprint(message)

            current_queue = self.model_retraining_queue.get()
            match_items = ['model_type', 'model_version', 'retrain_location']
            # Check for a duplicate entry in the retraining queue
            for task_item in current_queue['queue']:
                if all(task_item[key] == message[key] for key in match_items):
                    break
            else:
                # If the entry is a new entry we can look at the details, to start we need to look into the
                # model information collection in the database
                print('Requested model is not in the retraining queue')
                model_type = message['model_type']
                db_request = {'auth': {'port': 27017,
                                       'user': 'USER',
                                       'password': 'PASS'},
                              'request': {'request_type': 'query_data',
                                          'search_term': ['model_type', model_type],
                                          'search_type': 'model_information'}}

                # Look through the request returns
                database_ip_address = 'XX'
                database_port_number = 0
                server_url = 'http://%s:%s' % (database_ip_address, database_port_number)
                print('Sending request: %s' % db_request)
                db_response = requests.post(server_url, json=db_request,
                                            headers={'content-type': 'application/authjson'},
                                            timeout=10)
                response_content = json.loads(db_response.content.decode('utf_8'))
                pprint(response_content)
                pprint(db_request)

            response_code = 200
            self.send_response(response_code)
            self.send_header('Content_type', 'application/json')
            self.end_headers()
            self.wfile.write('testing'.encode(encoding='utf_8'))

    return RetrainingHandler


try:
    current_time = datetime.datetime.now()
    server_testing = True
    if server_testing:
        port_number = 9999
        http_address = '127.0.0.1'
    else:
        port_number = 20123
        http_address = '127.0.0.1'

    base_retraining_queue = queue.Queue()
    current_retraining_status = {'queue': [], 'running': []}
    #current_retraining_status['queue'].append({'model_type': 'photodeg_rf',
    #                                           'model_version': 'round5',
    #                                           'retrain_location': 'local',
    #                                           'new_dataset': {}})
    base_retraining_queue.put(current_retraining_status)
    retraining_manager = threading.Thread(target=retraining_thread,
                                          args=[base_retraining_queue])
    # retraining_manager.start()

    print('Starting the httpserver at %s' % current_time.strftime("%Y%m%d %H:%M:%S"))
    HandlerClass = MakeHandlerClass(base_retraining_queue)
    httpd = HTTPServer((http_address, port_number), HandlerClass)
    print('Started httpserver on port "%s" and httpaddress "%s"' % (port_number, http_address))
    httpd.serve_forever()

except KeyboardInterrupt:
    print('^C received, shutting down the web server')
    current_time = datetime.datetime.now()
    print('Stopping the httpserver at %s' % current_time.strftime("%Y%m%d %H:%M:%S"))
    httpd.socket.close()
