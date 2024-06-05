from intervaltree import Interval, IntervalTree
from typing import Type, TypeVar, List
from collections import OrderedDict


gtf = TypeVar("gtf", bound="GtfGff")


class GtfGff:

    def __init__(self):
        self._record_hashes = []
        self.records = OrderedDict()
        self.feature_index = {}
        self.seqname_index = {}
        self.attribute_index = {}
        self.metadata = {}

    def add_record(
        self,
        record: dict,
        linetype: str = "record",
        record_hash: int = None
    ):

        if linetype == "meta":
            self.metadata.update(record)
            return

        if not isinstance(record, dict):
            raise TypeError
        if record_hash is None:
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
        records = self._get_records(
            self._lookup_hash(self.feature_index, feature_type)
        )
        return records

    def get_records_by_seqname(self, seqname: str | int | list) -> list:

        if isinstance(seqname, int):
            seqname = str(seqname)
        records = self._get_records(
            self._lookup_hash(self.seqname_index, seqname)
        )
        return records

    def get_records_by_attribute(self, lookup_dict: dict) -> list:
        hashes = self._lookup_attribute_hashes(lookup_dict)
        return self._get_records(set(hashes))

    def __getitem__(self, idx: int | slice | list):

        if type(idx) not in (int, slice, list):
            raise TypeError(
                "Expected types int, slice, or list;"
                f"got '{idx}' of type: {type(idx)}"
            )

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

    def _process_query_and(self, value, depth=0):
        processed = self._process_query(value, depth=depth)
        processed = [set(lis) for lis in processed if isinstance(lis, list)]
        new_hashes = set(processed.pop(0))
        new_hashes = new_hashes.intersection(*processed)
        return new_hashes

    def _process_query_or(self, value, depth=0):
        processed = [self._process_query(v, depth=depth)[0] for v in value]
        return [item for sublist in processed for item in sublist]

    def _process_query_not(self, value, hashes, depth=0):

        processed = self._process_query(value, depth=depth)[0]
        if hashes == []:
            hashes = [self._record_hashes]
        new_hashes = []
        for sub in hashes:
            new_hashes.append(set(sub) - set(processed))

        return new_hashes

    def _process_query(self, conditions: dict, depth=0) -> list:

        depth += 1
        hashes = []
        try:
            for key, value in conditions.items():

                if key == "AND":
                    hashes.append(self._process_query_and(value, depth))

                elif key == "OR":
                    hashes.append(self._process_query_or(value, depth))

                elif key == "NOT":
                    hashes = self._process_query_not(value, hashes, depth)

                elif key == "seqname":
                    hashes.append(self._lookup_hash(self.seqname_index, value))

                elif key == "feature":
                    hashes.append(self._lookup_hash(self.feature_index, value))

                elif key == "attributes":
                    hashes.append(self._lookup_attribute_hashes(value))

                else:
                    raise ValueError(f"Invalid condition key: {key}")

            if depth == 1:
                hashes = hashes[0]

        except KeyError as e:
            print(
                f"KeyError: {e}."
                "Please check if the provided keys exist in the indices."
                )
            return []

        except Exception as e:
            print(f"An error occurred: {e}")
            return []

        return hashes

    def _process_query_string(self):
        pass

    def query(self, query: str | dict, return_records=False) -> gtf | dict:

        if isinstance(query, str):
            query = self._process_query_string(query)

        hashes = self._process_query(query)
        records = self._get_records(hashes)

        if return_records:
            return records
        else:
            return self.gtf_gff_from_records(records)

    def export_records(self):
        # list comprehension ensures a copy of each record is made
        return [dict(record) for record in self.records.values()]

    def _remove_record(self, record_hash):

        record = self.records[record_hash]
        self.feature_index[record["feature"]].remove(record_hash)
        self.seqname_index[record["seqname"]].remove(record_hash)
        for key, value in record["attributes"].items():
            self.attribute_index[key][value].remove(record_hash)
        self._record_hashes.remove(record_hash)
        del self.records[record_hash]

    def remove_empty_field(self, field: str | tuple | list):
        empty = []
        for h in self._record_hashes:
            try:
                if isinstance(field, str):
                    self.records[h][field]
                elif isinstance(field, tuple or list):
                    self.records[h][field[0]][field[0]]
            except KeyError:
                empty.append(h)

        for e in empty:
            self._remove_record(e)

    def __get_state__(self):
        return self.__dict__

    def __set_state(self, state):
        return self.__dict__.update(state)


def parse_gtf(filename, gtf=None):

    # if gtf is not give make one, if wrong type throw error
    if gtf is None:
        gtf = GtfGff()
    elif not isinstance(gtf, GtfGff):
        raise TypeError(f"gtf provided is not {GtfGff}, {type(gtf)} found\n")

    with open(filename, 'r') as file:
        for line in file:

            # Save comment lines as attributes
            if line.startswith('#'):
                record = line.lstrip("#!").rstrip("\n").split(" ")
                record = {record[0]: record[1]}
                gtf.add_record(record, linetype="meta")
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
            gtf.add_record(record, linetype="record")

    return gtf


def records_to_interval_tree(records: List[dict]) -> IntervalTree:
    interval_tree = IntervalTree()

    def add_inter(r):
        interval_tree.add(Interval(r.pop("start"), r.pop("end"), data=r))
    for r in records:
        add_inter(r)

    return interval_tree
