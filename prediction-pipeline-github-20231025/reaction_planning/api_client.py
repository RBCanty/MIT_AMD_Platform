import json
import requests
import time
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def parse_response(response):
    """Check response status and get JSON result."""
    try:
        result = response.json()
    except json.decoder.JSONDecodeError:
        result = None

    if response.status_code != 200:
        message = 'Request returned status code {0}.'.format(response.status_code)
        if result is not None and 'detail' in result:
            message += ' {0}'.format(result['detail'])
        print(message)

    return result


def verify_token(fn):
    """Check if the current token is expired or close to expiry."""

    def wrapper(self, *args, **kwargs):
        if self.token is not None:
            if time.time() > self.token_expiry:
                print('Your token has expired, please re-authenticate.')
            elif self.token_expiry - time.time() <= 60:
                # Less than a minute until expiration, automatically refresh
                self.refresh_token()
        return fn(self, *args, **kwargs)

    return wrapper


class APIClient(object):
    """
    Example of a simple API client for accessing the ASKCOS API.

    Usage:
        from api_client import APIClient
        client = APIClient('https://myhost.com')

        # list buyables
        result = client.get('buyables')

        # search for benzene in buyables
        result = client.get('buyables', params={'q': 'c1ccccc1'})

        # initiate retro task
        result = client.post('retro', data={'target': 'Brc1ccccc1'})

        # request token for authentication required endpoints
        client.authenticate('username', 'password')
    """

    def __init__(self, hostname, version=2, verify=True):
        self.client = requests.Session()
        self.client.verify = verify
        self.hostname = hostname
        self.version = version
        self.base_url = '{0}/api/v{1}'.format(self.hostname, self.version)
        self.token = None
        self.token_expiry = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.client.close()

    def close(self):
        """Close the current HTTP session."""
        self.client.close()

    def get_url(self, endpoint):
        """Return properly formatted URL based on the requested API endpoint."""
        return '{0}/{1}/'.format(self.base_url, endpoint.strip('/'))

    def update_auth_header(self, token):
        """Update authorization field in request header with token."""
        self.token = token
        self.client.headers.update({'Authorization': 'Bearer {}'.format(self.token)})

    def clear_token(self):
        """Clear authentication token and remove from request header."""
        self.token = None
        self.token_expiry = None
        del self.client.headers['Authorization']

    def refresh_token(self):
        """Refresh an unexpired token by obtaining a new token."""
        self.token_expiry = time.time() + 300  # Get time first to be conservative
        result = self.post('/token-refresh/', data={'token': self.token})
        self.update_auth_header(result['token'])

    def authenticate(self, username, password):
        """Request a new token using username and password."""
        self.token_expiry = time.time() + 300  # Get time first to be conservative
        result = self.post('/token-auth/', data={'username': username, 'password': password})
        self.update_auth_header(result['token'])

    @verify_token
    def get(self, endpoint, params=None):
        """Send a get request with optional query parameters."""
        resp = self.client.get(self.get_url(endpoint), params=params)
        return parse_response(resp)

    @verify_token
    def post(self, endpoint, data=None):
        """Send a post request with optional data payload."""
        resp = self.client.post(self.get_url(endpoint), data=data)
        return parse_response(resp)

    @verify_token
    def options(self, endpoint):
        """Send an options request to get parameter information."""
        resp = self.client.options(self.get_url(endpoint))
        return parse_response(resp)

    def get_result(self, task_id, timeout=600, interval=5):
        """Retrieve celery task output"""
        for _ in range(timeout // interval):
            result = self.get('/celery/task/{0}/'.format(task_id))
            if result.get('complete') or result.get('failed'):
                return result
            else:
                time.sleep(interval)

    def batch_get(self, batch):
        """
        Send a batch of get requests.

        Accepts a list of (endpoint, params) tuples.
        """
        return [self.get(endpoint, params=params) for endpoint, params in batch]

    def batch_post(self, batch):
        """
        Send a batch of post requests.

        Accepts a list of (endpoint, data) tuples.
        """
        return [self.post(endpoint, data=data) for endpoint, data in batch]

    def batch_get_result(self, task_ids, timeout=600, interval=5):
        """Retrieve celery task output for a batch of tasks."""
        return [self.get_result(task_id, timeout=timeout, interval=interval) for task_id in task_ids]


if __name__ == '__main__':

    client = APIClient('https://localhost', verify=False)
    client.authenticate('testuser', 'reallybadpassword')
    result = client.get('results')
    print(result)
