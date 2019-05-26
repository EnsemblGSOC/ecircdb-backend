# ecircdb-backend
Backend for the ecircdb

## Setup on your local machine

1. Clone the Repository
```shell
$ git clone https://github.com/EnsemblGSOC/ecircdb-backend.git ecircdb-backend
$ cd ecircdb-backend
```
2. Make virtual environment setup
```shell
$ pip install virtualenv
$ virtualenv -p python3 <env_name>
```
3. Activate your environment
```shell
$ source <env_name>/bin/activate
```
4. Install requirements
```shell
$ pip install -r requirement.txt
```
5. Migrate Files
```shell
$ python manage.py makemigrations
$ python manage.py migrate
```
6. Create a superuser
```shell
$ python manage.py createsuperuser
```
6. Start the server
```shell
$ python manage.py runserver
```
7. Visit `localhost:8000/control` and login using the account created in above step. This is the admin panel from where data can be added, deleted or modified.
8. Check other urls from the list at `localhost:8000`