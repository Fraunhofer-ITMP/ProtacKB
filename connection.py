# -*- coding: utf-8 -*-

import logging

from py2neo import Graph, SystemGraph
#from py2neo.database.work import ClientError

from constants import FRAUNHOFER_ADMIN_NAME, FRAUNHOFER_ADMIN_PASS, FRAUNHOFER_URL

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)

system_graph = SystemGraph(FRAUNHOFER_URL, auth=(FRAUNHOFER_ADMIN_NAME, FRAUNHOFER_ADMIN_PASS))


# def create_new_user():
#     """Creates a new user with public roles in the DB."""
#
#     print("SUCCESS: Connected to the Neo4j Database.")
#
#     # TODO: Get name list from the survey sheet
#
#     for name in ["test"]:
#         try:
#             system_graph.run(
#                 """CALL dbms.security.createUser("{}", "pass");""".format(name)
#             )
#             logger.info(f"Added new user {name}")
#         except ClientError:
#             logger.info("User exists!")


def create_new_db(db_name: str):
    """Create a new DB table"""

    system_graph.run("""CREATE OR REPLACE DATABASE {}""".format(db_name))

    logger.info(f"Created {db_name} in Neo4J!")


def check_database(db_name: str) -> bool:
    database_info = system_graph.run("""SHOW DATABASES""").data()

    for info in database_info:
        if info["name"] == db_name:
            return True
    return False


def populate_db(db_name: str):
    db_name=db_name.lower()
    in_db = check_database(db_name)

    if not in_db:
        create_new_db(db_name)
        in_db = check_database(db_name)  # to update the database

    assert in_db is True

    conn = Graph(URL, auth=(ADMIN_NAME, ADMIN_PASS), name=db_name)
    conn.delete_all()
    tx = conn.begin()
    return tx,conn


def _add_nodes(
    node_dict: dict,
    tx
):
    """Add nodes from dictionary."""
    for node_type in node_dict:
        tx.create(node_dict[node_type])