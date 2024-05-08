import re
from typing import Type, TypeVar
from collections import OrderedDict

gtf = TypeVar("gtf", bound = "GtfGff")

class GtfGff:


    def __init__(self):
        self._record_hashes = []
        self.records = OrderedDict()
        self.feature_index = {}
        self.seqname_index = {}
        self.attribute_index = {}
        self.metadata = {}
    
    def add_record(self, record: dict, linetype: str = "record", record_hash: int = None):

        if linetype == "meta":
            self.metadata.update(record)
            return

        if type(record) != dict:
            raise TypeError
        if  record_hash == None:
            record_hash = hash(str(record))
        
        self._record_hashes.append(record_hash)
        self.records[record_hash] = record
            
        feature_type = record["feature"]
        if feature_type not in self.feature_index.keys():
            self.feature_index[feature_type] = []
        self.feature_index[feature_type].append(record_hash)

        seqname = record["seqname"]
        if seqname not in self.seqname_index.keys():
            self.seqname_index[seqname] = []
        self.seqname_index[seqname].append(record_hash)

        # Indexing by attributes
        for attribute, value in record['attributes'].items():
            if attribute not in self.attribute_index:
                self.attribute_index[attribute] = {}
            if value not in self.attribute_index[attribute]:
                self.attribute_index[attribute][value] = []
            self.attribute_index[attribute][value].append(record_hash)
    
    def get_records_by_feature(self, feature_type):
        hashes = self.feature_index.get(feature_type, [])
        return [self.records[h] for h in hashes]
    
    @staticmethod
    def _lookup_hash(index: dict, keys: str | list):
        if keys is None:
            return []
        if isinstance(keys, str):
            return index.get(keys, [])
        else:
            return list(set(h for key in keys for h in index.get(key)))
    
    def _lookup_attribute_hashes(self, lookup_dict):
        hashes = []
        for k, v in lookup_dict.items():
            hashes += self._lookup_hash(self.attribute_index[k], v)
        return hashes

    
    def _get_records(self, hashes: int | list | set) -> list:
        if isinstance(hashes, int):
            return [self.records[hashes]]
        else:
            return [self.records[h] for h in hashes]
    

    def get_records_by_feature(self, feature_type: str | list) -> list:
        return self._get_records(self._lookup_hash(self.feature_index, feature_type))

    def get_records_by_seqname(self, seqname: str | list) -> list:
        return self._get_records(self._lookup_hash(self.seqname_index, seqname))

    def get_records_by_attribute(self, lookup_dict: dict) -> list:
        hashes = self._lookup_attribute_hashes(lookup_dict)
        return self._get_records(set(hashes))

        
    def __getitem__(self, idx: int | slice | list ):

        if type(idx) not in (int, slice, list):
            raise TypeError(f"Expected types int, slice, or list; got '{idx}' of type: {type(idx)}") 

        if isinstance(idx, int):
            return self.records[self._record_hashes[idx]]

        if isinstance(idx, slice):
            records = [
                self.records[self._record_hashes[i]] for i in range(
                    idx.start,
                    idx.stop,
                )
            ] 
            return records
        
        if isinstance(idx, list):
            return [self.records[self._record_hashes[i]] for i in idx]

    def __len__(self):
        return len(self._record_hashes)
    
    @classmethod
    def gtf_gff_from_records(cls: Type[gtf], records) -> gtf:

        new_gtf = cls()
        if isinstance(records, dict):
            new_gtf.add_record(records)
        else:
            for r in records:
                new_gtf.add_record(r)
        return new_gtf
    

    def _process_query(self, conditions: dict, count = 0) -> list:

        count += 1
        hashes = []
        for key, value in conditions.items():

            if key == "AND":
                processed = self._process_query(value, count = count)
                new_hashes = set(processed.pop(0))
                hashes.append(new_hashes.intersection(*processed))

            elif key == "OR":
                processed = [self._process_query(v, count = count)[0] for v in value] 
                processed = [item for sublist in processed for item in sublist]
                hashes.append(processed)

            elif key == "NOT":
                processed = self._process_query(value, count = count)[0]
                new_hashes = []
                for sub in hashes:
                    new_hashes.append(set(sub) - set(processed))
                hashes = new_hashes

            if key == "seqname":
                hashes.append(self._lookup_hash(self.seqname_index, value))
            elif key == "feature":
                hashes.append(self._lookup_hash(self.feature_index, value))
            elif key == "attributes":
                hashes.append(self._lookup_attribute_hashes(value))

        if count == 1:
            hashes = hashes[0]
        
        return hashes



    
def parse_gtf(filename, gtf = None):

    """
    Read GTF (Gene Transfer Format) or GFF (General Feature Format) files and parse the content.

    Parameters:
        file_path (str): The path to the GTF/GFF file.

    Returns:
        GtfGff data structure defined above
        
    """
    # if gtf is not give make one, if wrong type throw error
    if gtf == None:
        gtf = GtfGff()
    else:
        assert type(gtf) == GtfGff, f"gtf provided is not {GtfGff}, {type(gtf)} found\n"
    
    with open(filename, 'r') as file:
        for line in file:

            # Save comment lines as attributes
            if line.startswith('#'):
                record = line.lstrip("#!").rstrip("\n").split(" ")
                record = {record[0]: record[1]}
                gtf.add_record(record, linetype = "meta")
                continue

            # Split line by tab
            fields = line.strip().split('\t')
            # Parse attributes
            attributes = {}
            for attribute in fields[8].split(';'):
                key_value = attribute.strip().split(' ')
                if len(key_value) == 2:
                    attributes[key_value[0]] = key_value[1].strip('"')

            # Create feature dictionary
            record = {
                'seqname': fields[0],
                'source': fields[1],
                'feature': fields[2],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'score': fields[5],
                'strand': fields[6],
                'frame': fields[7],
                'attributes': attributes
            }
            gtf.add_record(record, linetype = "record")
    
    return gtf

