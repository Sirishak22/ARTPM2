##############################################
# Filter CUIs in termfile by number of atoms must be greater than 1. Requires local UMLS database.
# Lowell Milliken
#
#
##############################################
import sys
import pymysql.cursors
import pickle

def cui_filter(termfile):
    host = 'localhost'
    user = 'umlsuser'
    passw = 'umlspass'
    db = 'umls'
    charset = 'utf8'

    connection = pymysql.connect(host=host,
                                 user=user,
                                 password=passw,
                                 db=db,
                                 charset=charset,
                                 cursorclass=pymysql.cursors.DictCursor)

    with open(termfile, 'rb') as infile:
        terms = pickle.load(infile)

    fterms = []

    try:
        sql = 'SELECT str FROM mrconso' \
              ' WHERE cui = %s' \
              ' AND ispref = \'Y\'' \
              ' AND ts = \'P\'' \
              ' AND stt = \'PF\'' \
              ' AND lat = \'ENG\';'

        for cui in terms:
            with connection.cursor() as cursor:
                cursor.execute(sql, (cui,))
                result = cursor.fetchone()
                if result is not None and len(result['str'].split()) > 1:
                    print(result)
                    fterms.append(cui)
    finally:
        connection.close()

    return fterms

if __name__ == '__main__':
    cui_filter(sys.argv[1])