# from sqlalchemy import create_engine, Column, Integer, String, Boolean, ForeignKey

import datetime
import os
from abc import abstractmethod
from typing import List, Optional

from sqlalchemy import (Boolean, Column, Integer, MetaData, String, Table,
                        create_engine, text)


class SoftwareItem:
    def __init__(
        self, name, path, database, installed, env_path, 
        tag: str = "undefined", db_version: Optional[str] = None, 
        needs_update: bool = False, binary_name: Optional[str] = None
    ) -> None:
        self.name = name
        self.path = path
        self.database = database
        self.installed = installed
        self.env_path = env_path
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.tag = tag
        self.db_version = db_version
        self.needs_update = needs_update
        self.binary_name = binary_name

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.database}, {self.installed}, {self.env_path})"


class DatabaseItem:
    def __init__(self, name, path, installed, software: str = "none",
                 version: Optional[str] = None, source_url: Optional[str] = None, 
                 file_mod_date: Optional[str] = None, description: Optional[str] = None,
                 db_category: Optional[str] = None, db_name: Optional[str] = None,
                 db_type: Optional[str] = None) -> None:
        self.name = name
        self.path = path
        self.installed = installed
        self.software = software
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.version = version
        self.source_url = source_url
        self.file_mod_date = file_mod_date
        self.description = description
        self._parse_name_fields(name, db_category, db_name, db_type)

    def _parse_name_fields(self, name: str, db_category: Optional[str], 
                           db_name: Optional[str], db_type: Optional[str]):
        if '/' in name and db_category is None:
            parts = name.split('/', 1)
            self.db_category = parts[0]
            self.db_name = parts[1]
        else:
            self.db_name = db_name if db_name else name
            self.db_category = db_category if db_category else self.software
        self.db_type = db_type if db_type else self.software

    def __repr__(self) -> str:
        return f"({self.name}, {self.path}, {self.installed})"



class Utility_Repository:
    """Communicates with sql database to add a software dbs"""

    database_item = SoftwareItem
    software_item = DatabaseItem
    dbtype_local: str = "sqlite"
    SOFTWARE_TABLE_NAME: str = "software"
    DATABASE_TABLE_NAME: str = "database"
    
    tables: list = ["software", "database"]

    def __init__(self, db_path="", install_type="local") -> None:
        self.db_path = db_path

        self.setup_engine(install_type)

        # self.connection = self.engine.connect()
        self.metadata = MetaData()

        if not self.check_tables_exists():
            # self.delete_tables()
            # self.clear_tables()
            self.create_tables()

    def setup_engine(self, install_type):
        if not os.path.exists(self.db_path):
            os.makedirs(self.db_path, exist_ok=True)
        if install_type == "local":
            self.setup_engine_local()
        elif install_type == "docker":
            self.setup_engine_docker()

    def setup_engine_local(self):
        self.engine = create_engine(
            f"{self.dbtype_local}:////"
            + os.path.join(*self.db_path.split("/"), "utility_local.db")
        )

    def setup_engine_docker(self):
        self.engine = create_engine(
            f"{self.dbtype_local}:////"
            + os.path.join(*self.db_path.split("/"), "utility_local.db")
        )

    def setup_engine_posrgres(self):
        from decouple import config

        self.engine = create_engine(
            f"postgresql+psycopg2://{config('DB_USER')}:{config('DB_PASSWORD')}@{config('DB_HOST')}:{config('DB_PORT')}/{config('DB_NAME')}"
        )

    def engine_execute_return_table(self, string: str):
        sql = text(string)

        rows = None

        with self.engine.connect() as conn:
            result = conn.execute(sql)
            #conn.commit()
            rows = result.fetchall()

        return rows

    def check_table_exists(self, table_name):
        """
        Check if a table exists in the database
        """
        find = self.engine.execute(
            f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}'"
        ).fetchall()

        if len(find) > 0:
            return True
        else:
            return False

    def check_tables_exists(self):
        """
        Check if the tables exist in the database
        """
        for table_name in self.tables:
            if not self.check_table_exists(table_name):
                return False
        return True

    def create_software_table(self):
        self.software = Table(
            self.SOFTWARE_TABLE_NAME,
            self.metadata,
            Column("name", String, primary_key=True),
            Column("path", String),
            Column("database", String),
            Column("installed", Boolean),
            Column("tag", String, default="undefined"),
            Column("env_path", String),
            Column("date", String),
            Column("db_version", String),
            Column("needs_update", Boolean),
        )

        self.engine.execute(
            "CREATE TABLE IF NOT EXISTS software (name TEXT PRIMARY KEY, path TEXT, database TEXT, installed BOOLEAN, tag TEXT, env_path TEXT, date TEXT, db_version TEXT, needs_update BOOLEAN)"
        )

    def create_database_table(self):
        self.database = Table(
            self.DATABASE_TABLE_NAME,
            self.metadata,
            Column("id", Integer, primary_key=True, autoincrement=True),
            Column("name", String),
            Column("db_category", String),
            Column("db_name", String),
            Column("db_type", String),
            Column("path", String),
            Column("installed", Boolean),
            Column("software", String),
            Column("date", String),
            Column("version", String),
            Column("source_url", String),
            Column("file_mod_date", String),
            Column("description", String),
        )

        self.engine.execute(
            "CREATE TABLE IF NOT EXISTS database (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, db_category TEXT, db_name TEXT, db_type TEXT, path TEXT, installed BOOLEAN, software TEXT, date TEXT, version TEXT, source_url TEXT, file_mod_date TEXT, description TEXT)"
        )

    def delete_tables(self):
        self.delete_table("software")
        self.delete_table("database")

    def delete_table(self, table_name):
        self.engine.execute(f"DROP TABLE {table_name}")

    def clear_tables(self):
        self.clear_table("software")
        self.clear_table("database")

    def clear_table(self, table_name):
        self.engine.execute(f"DELETE FROM {table_name}")

    def print_table_schema(self, table_name):
        print(self.engine.execute(f"PRAGMA table_info({table_name})").fetchall())

    def create_tables(self):
        """
        Create the tables
        """

        self.create_software_table()
        self.create_database_table()

        self.metadata.create_all(self.engine)

    def dump_software(self, directory: str):
        """
        Dump the software table to a tsv file
        """

        self.dump_table_tsv("software", directory)

    def dump_database(self, directory: str):
        """
        Dump the database table to a tsv file
        """

        self.dump_table_tsv("database", directory)

    def dump_tables(self, directory: str):
        """
        Dump the database & software tables to a tsv file
        """

        self.dump_table_tsv("software", directory)
        self.dump_table_tsv("database", directory)

    def dump_table_tsv(self, table_name: str, directory: str):
        """
        Dump a table to a tsv file
        """

        if table_name not in self.tables:
            print(f"Table {table_name} not found. Available tables: {self.tables}")
            return

        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)

        with open(os.path.join(directory, f"{table_name}.tsv"), "w") as f:
            for row in self.engine.execute(f"SELECT * FROM {table_name}"):
                f.write("\t".join([str(x) for x in row]) + "\n")

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

    def select_explicit_statement(
        self, table_name, field, id, filters: List[tuple] = []
    ):
        """
        select from table.
        """
        sql_statement = f"SELECT * FROM {table_name} WHERE {field}='{id}'"
        for filter in filters:
            column_name, value = filter

            if self.check_column_exists(table_name, column_name) is False:
                continue
            if value is None:
                continue

            sql_statement += f" AND {column_name}='{value}'"

        return sql_statement

    def check_column_exists(self, table_name, column_name):
        """
        Check if a column exists in a table
        """

        from sqlalchemy import inspect

        inspector = inspect(self.engine)
        columns = inspector.get_columns(table_name)
        find = any([i["name"] == column_name for i in columns])

        if find:
            return True
        else:
            return False

    def check_exists(self, table_name: str, id: str):
        """
        Check if a record exists in a table
        """

        check_list = [id]
        if "_" in id:
            check_list.append(id.split("_")[0])
        check_list = [f"'{i}'" for i in check_list]
        check_list = ",".join(check_list)

        find = self.engine_execute_return_table(
            f"SELECT * FROM {table_name} WHERE name='{id}'"
        )

        find = len(find) > 0
        if find:
            return True
        else:
            return False

    @abstractmethod
    def add_software(self, item: software_item):
        """
        Add a record to a table
        """

        self.engine.execute(
            f"INSERT INTO software (name, path, database, installed, env_path, date) VALUES ('{item.name}', '{item.path}', '{item.database}', '{item.installed}', '{item.env_path}', '{item.date}')"
        )

    @abstractmethod
    def add_database(self, item: database_item):
        """
        Add a record to a table
        """

        self.engine.execute(
            f"INSERT INTO database (name, path, installed, date) VALUES ('{item.name}', '{item.path}', '{item.installed}', '{item.date}')"
        )
