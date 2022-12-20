# from sqlalchemy import create_engine, Column, Integer, String, Boolean, ForeignKey

import os
from abc import abstractmethod

from sqlalchemy import Boolean, Column, MetaData, String, Table, create_engine


class software_item:
    def __init__(self, name, path, database, installed, env_path) -> None:
        self.name = name
        self.path = path
        self.database = database
        self.installed = installed
        self.env_path = env_path

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.database}, {self.installed}, {self.env_path})"


class database_item:
    def __init__(self, name, path, installed, software: str = "none") -> None:
        self.name = name
        self.path = path
        self.installed = installed
        self.software = software

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.installed})"


class Utility_Repository:
    """Communicates with sql database to add a software dbs"""

    database_item = database_item
    software_item = software_item
    dbtype_local: str = "sqlite"

    def __init__(self, db_path="", install_type="local") -> None:
        self.db_path = db_path

        self.setup_engine(install_type)

        # self.connection = self.engine.connect()

        self.metadata = MetaData()
        self.create_tables()

    def setup_engine(self, install_type):
        if install_type == "local":
            self.setup_engine_local()
        elif install_type == "docker":
            self.setup_engine_docker()

    def setup_engine_local(self):
        self.engine = create_engine(
            f"{self.dbtype_local}:////"
            + os.path.join(*self.db_path.split("/"), "utility.db")
        )

    def setup_engine_docker(self):

        self.engine = create_engine(
            f"{self.dbtype_local}:////"
            + os.path.join(*self.db_path.split("/"), "utility.db")
        )

    def setup_engine_posrgres(self):
        from decouple import config

        self.engine = create_engine(
            f"postgresql+psycopg2://{config('DB_USER')}:{config('DB_PASSWORD')}@{config('DB_HOST')}:{config('DB_PORT')}/{config('DB_NAME')}"
        )

    def create_software_table(self):

        self.software = Table(
            "software",
            self.metadata,
            Column("name", String),
            Column("path", String),
            Column("database", String),
            Column("installed", Boolean),
            Column("env_path", String),
        )

    def create_database_table(self):
        self.database = Table(
            "database",
            self.metadata,
            Column("name", String),
            Column("path", String),
            Column("installed", Boolean),
            Column("software", String),
        )

    def create_tables(self):
        self.create_software_table()
        self.create_database_table()

        self.metadata.create_all(self.engine)

    def get_by_name(self, table_name, id):
        """
        Get a record by id from a table
        """

        return self.engine.execute(f"SELECT * FROM {table_name} WHERE name='{id}'")

    def select_explicit(self, table_name, field, id):
        """
        select from table.
        """
        sql_statement = f"SELECT * FROM {table_name} WHERE {field}='{id}'"

        find = self.engine.execute(sql_statement)

        return find

    def get_list_tables(self):
        """
        Get a list of tables
        """

        find = self.engine.execute("SELECT name FROM sqlite_master WHERE type='table'")
        find = [i[0] for i in find]
        return find

    def get_list_unique_field(self, table_name, id):
        """
        Get a list of unique values in a field
        """

        find = self.engine.execute(f"SELECT DISTINCT {id} FROM {table_name}")

        find = [i[0] for i in find]
        return find

    def select_explicit_statement(self, table_name, field, id):
        """
        select from table.
        """
        sql_statement = f"SELECT * FROM {table_name} WHERE {field}='{id}'"

        return sql_statement

    def check_exists(self, table_name, field, id):
        """
        Check if a record exists in a table
        """

        check_list = [id]
        if "_" in id:
            check_list.append(id.split("_")[0])
        check_list = [f"'{i}'" for i in check_list]
        check_list = ",".join(check_list)

        find = self.engine.execute(
            f"SELECT * FROM {table_name} WHERE {field} IN ({check_list})"
        ).fetchall()
        find = len(find) > 0

        if find:
            return True
        else:
            return False

    def add_software(self, item):
        """
        Add a record to a table
        """

        self.engine.execute(
            f"INSERT INTO software (name, path, database, installed, env_path) VALUES ('{item.name}', '{item.path}', '{item.database}', '{item.installed}', '{item.env_path}')"
        )

    @abstractmethod
    def add_database(self, item):
        """
        Add a record to a table
        """

        self.engine.execute(
            f"INSERT INTO database (name, path, installed) VALUES ('{item.name}', '{item.path}', '{item.installed}')"
        )
