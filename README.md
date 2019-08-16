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

4. Install `libmysqlclient`

For example on Debian/Ubuntu you must install the package:

```shell
sudo apt-get install libmysqlclient-dev
```

For recent versions of debian/ubuntu (as of 2018) it is:

```shell
sudo apt install default-libmysqlclient-dev
```

5. Install requirements

```shell
$ pip install -r requirement.txt
```

6. Create a MySQL user and database normally, remember the name of the database and credentials for the user
7. Configure the .env file for the project, you can clone the available `.en.stencil` file to checkout the format for .env file

```shell
$ cp .env.stencil .env
```

8. Fill out the database name, user, password in the .env file, you can specify the host and port for the database optionally
9. Migrate Files

```shell
$ python manage.py makemigrations
$ python manage.py migrate
```

10. Create a superuser

```shell
$ python manage.py createsuperuser
```

11. Start the server

```shell
$ python manage.py runserver
```

12. Visit `localhost:8000/control` and login using the account created in above step. This is the admin panel from where data can be added, deleted or modified.
13. Check other urls from the list at `localhost:8000`
