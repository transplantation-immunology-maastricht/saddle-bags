# This file is part of saddle-bags.
#
# saddle-bags is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# saddle-bags is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with saddle-bags. If not, see <http://www.gnu.org/licenses/>.

from saddlebags.AlleleSubCommon import resourcePath

import sqlite3
from sqlite3 import Error

from BioSQL import BioSeqDatabase

import logging

from os.path import join
from sys import exc_info

from Bio.SeqIO import parse

def connectSqliteDatabase(databaseFullPath):
    # Create a sqlite database
    conn = None
    try:
        conn = sqlite3.connect(databaseFullPath)
        logging.debug('Sqlite version=' + str(sqlite3.version))
        logging.info('Created connection to database:' + str(databaseFullPath))
        return conn
    except Error as e:
        logging.error(e)
        if conn:
            conn.close()
    #finally:
#        if conn:
#            conn.close()


# Silly Create Table method
def createTable(databaseConnection, tableName):
    sqlStatement = ('CREATE TABLE IF NOT EXISTS ' + tableName + ' ( '
    + 'id integer PRIMARY KEY,'
    + 'name text NOT NULL,'
    + 'begin_date text,'
    + 'end_date text'
    + '); ')

    logging.debug('Creating Table!')
    executeSqlStatement(databaseConnection, sqlStatement)

#def queryTable(databaseConnection, tableName):


def executeSqlStatement(databaseConnection, sqlStatement):
    logging.debug('Executing Sql Statement...')
    try:
        c = databaseConnection.cursor()
        c.execute(sqlStatement)
    except Error as e:
        logging.error(e)

def executeSqlScript(databaseConnection, sqlScript):
    # I can call executescript() or execute(). execute() only allows a single statement.
    # I think because of security reasons.
    logging.debug('Executing Sql Script...')
    try:
        c = databaseConnection.cursor()
        c.executescript(sqlScript)
    except Error as e:
        logging.error(e)


def setupBioSqlDatabase(databaseFullPath, bioSqlCommandFilePath):
    logging.debug('Initializing the Biosql database!')
    # Create the database file.
    sqliteConnection=connectSqliteDatabase(databaseFullPath)

    # Run the command to initialize BioSql database
    bioSqlCommandFile = open(bioSqlCommandFilePath, 'r')

    bioSqlText = bioSqlCommandFile.read()
    executeSqlScript(sqliteConnection, bioSqlText)


    #server = BioSeqDatabase.open_database(driver='sqlite')

    # Import the hla.dat file.

    sqliteConnection.close()


def loadHLADataIntoBioSql(databaseFullPath, hlaDataFolder):
    # This code is adopted from the script included with SeqAnn
    # https://github.com/nmdp-bioinformatics/seq-ann/scripts/create_imgtdb.py
    # Which is also released under the GPL 3.0 license.
    logging.debug('Loading HLA data from this folder:' + hlaDataFolder)



    #server = BioSeqDatabase.open_database(driver="pymysql", user="root",
    #                                          passwd="", host="localhost",
    #                                             db="bioseqdb")

    sqliteConnection = connectSqliteDatabase(databaseFullPath)

    #for dbv in dblist:

    #hladat = download_dat(dbv)
    #allele_list = download_allelelist(dbv)
    hladat = join(hlaDataFolder, 'hla.dat')
    allele_list = join(hlaDataFolder, 'Allelelist.txt')

    hla_names = {}
    try:
        #s = "," if dbv == "3260" or dbv == "3270" else " "
        s = ','
        with open(allele_list, 'r') as f:
            for line in f:
                line = line.rstrip()
                accession, name = line.split(s)
                hla_names.update({accession: name})
        f.close()
        logging.debug("Loaded allele names " +  allele_list)
    except ValueError as err:
        logging.error("Allelelist error: {0}".format(err))
        sqliteConnection.close()
        #os.remove(hladat)
        #os.remove(allele_list)
        #sys.exit()

    try:
        seq_list = parse(hladat, "imgt")
    except ValueError as err:
        logging.error("Read dat error: {0}".format(err))
        sqliteConnection.close()
        #os.remove(hladat)
        #os.remove(allele_list)
        #sys.exit()

    new_seqs = {"A": [], "B": [], "C": [], "DRB1": [],
                "DQB1": [], "DRB3": [], "DRB4": [], "DRB5": [],
                "DQA1": [], "DPA1": [], "DPB1": [], "DRA": []}

    for seq in seq_list:
        if seq.name in hla_names:
            loc, allele = hla_names[seq.name].split("*")
            if loc in new_seqs:
                hla_name = "HLA-" + hla_names[seq.name]
                seq.name = hla_name
                new_seqs[loc].append(seq)

    #dbsp = list(dbv)
    #descr = ".".join([dbsp[0], dbsp[1] + dbsp[2], dbsp[3]])
    #print("Loaded IMGT dat file ", descr)

    for locus in new_seqs:
        #dbname = dbv + "_" + locus
        dbname = locus
        dbdescription = "IMGT/HLA " + locus
        db = sqliteConnection.new_database(dbname, description=dbdescription)
        try:
            count = db.load(new_seqs[locus])
        except:
            logging.error("Failed to load :" + str(exc_info()[0]))
            sqliteConnection.close()
            #os.remove(hladat)
            #os.remove(allele_list)
            #sys.exit()

        print("Loaded ", count, " for ", locus)
        sqliteConnection.commit()

    #os.remove(hladat)
    #os.remove(allele_list)
    logging.debug("Finished ", dbdescription)

    #server.close()
    sqliteConnection.close()

