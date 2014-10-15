import os
from datetime import datetime
import json
import boto
from boto.s3.key import Key
from bson.objectid import ObjectId
from pymongo import MongoClient
from six.moves.urllib.request import urlopen
from molecule import qcheminpfile, protonate, cid2mol
from fabric.api import local, task

#################################################################

JOBLOCK = 'current-job.lock'
INPUTFILE = 'qchemfile.dat'
OUTPUTFILE = 'calculation.out'

PUBCHEM_ASSAY = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1996/cids/json'
BUCKET_NAME = 'rmcgibbo-bagofwaves'

#################################################################

with open(os.path.join(os.path.dirname(__file__), '.env')) as f:
    for line in f:
        line = line.strip()
        if line:
            key, value = line.split(' = ')
            os.environ[key] = value.strip()


def connect_mongo():
    uri = 'mongodb://%s:%s@%s' % (os.environ['MONGO_USER'],
                                  os.environ['MONGO_PASSWORD'],
                                  os.environ['MONGO_URL'])
    client = MongoClient(uri, safe=True)
    cursor = client.bagofwaves.experiment
    return cursor

@task
def interact():
    import IPython as ip
    ip.embed()

def upload_to_s3(bucket, prefix, filename):
    k = Key(bucket)
    if not filename.endswith('.gz'):
        local('gzip %s' % filename)
        filename = filename + '.gz'
    k.key = os.path.join(prefix, filename)
    k.set_contents_from_filename(filename)


#################################################################


@task
def initialize_database():
    records = json.load(urlopen(PUBCHEM_ASSAY))
    cids = records['InformationList']['Information'][0]['CID']
    cids = cids[0:10]

    cursor = connect_mongo()
    for i, cid in enumerate(cids):
        cursor.insert({'cid': cid, 'status': 'NEW'})

        # DEBUGGING
        if i >= 10:
            break

@task
def generate_qchemfile():
    def get_cid():
        cursor = connect_mongo()
        record = cursor.find_and_modify(
            query={"status": "NEW"},
            update={"$set": {"status": "PENDING",
                             "started": datetime.now()}})

        if record is None:
            exit(1)
        return record

    record = get_cid()
    cid = record['cid']
    qcfile = qcheminpfile(protonate(cid2mol(cid)), comment=cid)

    with open(JOBLOCK, 'w') as f:
        f.write(str(record['_id']))

    with open(INPUTFILE, 'w') as f:
        f.write(qcfile)


@task
def run_qchem():
    local('qchem qchemfile.dat calculation.out')


@task
def upload():
    with open(JOBLOCK) as f:
        jobid = f.read()

    conn = boto.connect_s3(os.environ['AWS_ACCESS_KEY_ID'],
                           os.environ['AWS_SECRET_ACCESS_KEY'])
    bucket = conn.get_bucket(BUCKET_NAME)

    upload_to_s3(bucket, jobid, OUTPUTFILE)
    upload_to_s3(bucket, jobid, INPUTFILE)

    cursor = connect_mongo()
    cursor.find_and_modify(
        query={"_id": ObjectId(jobid)},
        update={"$set": {"status": "COMPLETE",
                         "completed": datetime.now()}})
