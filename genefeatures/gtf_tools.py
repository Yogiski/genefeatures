from __future__ import annotations
import re
from intervaltree import Interval, IntervalTree
from typing import List, Dict
from collections import OrderedDict


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
            raise TypeError(
                f"record argument must be type {dict}; "
                f"got {record} of type {type(record)}"
            )
        if record_hash is None:
            record_hash = hash(str(record))

        self._record_hashes.append(record_hash)
        self.records[record_hash] = record

        feature_type = record["feature"]
        self.feature_index.setdefault(feature_type, []).append(record_hash)

        seqname = record["seqname"]
        self.seqname_index.setdefault(seqname, []).append(record_hash)

        # Indexing by attributes
        for attribute, value in record['attributes'].items():
            self.attribute_index \
                .setdefault(attribute, {}) \
                .setdefault(value, []) \
                .append(record_hash)

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

    def _get_records(self, hashes: int | list | set) -> List[dict]:
        if isinstance(hashes, int):
            return [self.records[hashes]]
        else:
            return [self.records[h] for h in hashes]

    def get_records_by_feature(self, feature_type: str | list) -> List[dict]:
        records = self._get_records(
            self._lookup_hash(self.feature_index, feature_type)
        )
        return records

    def get_records_by_seqname(self, seqname: str | int | list) -> List[dict]:

        if isinstance(seqname, int):
            seqname = str(seqname)
        records = self._get_records(
            self._lookup_hash(self.seqname_index, seqname)
        )
        return records

    def get_records_by_attribute(self, lookup_dict: dict) -> List[dict]:
        hashes = self._lookup_attribute_hashes(lookup_dict)
        return self._get_records(set(hashes))

    def __getitem__(self, idx: int | slice | list) -> dict | List[dict]:

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
    def gtf_gff_from_records(cls: GtfGff, records: dict) -> GtfGff:

        new_gtf = cls()
        if isinstance(records, dict):
            new_gtf.add_record(records)
        elif isinstance(records, list):
            for r in records:
                new_gtf.add_record(r)
        else:
            raise TypeError(
                f"records must be type {dict}; "
                f"got {records} of type {type(records)}"
            )

        return new_gtf

    def _process_query_and(self, value, depth: int = 0):
        processed = self._process_query(value, depth=depth)
        processed = [set(lis) for lis in processed if isinstance(lis, list)]
        new_hashes = set(processed.pop(0))
        new_hashes = new_hashes.intersection(*processed)
        return new_hashes

    def _process_query_or(self, value, depth=0):
        processed = [self._process_query(v, depth=depth)[0] for v in value]
        return [item for sublist in processed for item in sublist]

    def _process_query_not(
        self,
        value: Dict[str, str],
        hashes: List[str],
        depth: int = 0
    ) -> List[str]:

        processed = self._process_query(value, depth=depth)[0]
        if not hashes:
            hashes = [self._record_hashes]
        new_hashes = [set(sub) - set(processed) for sub in hashes]

        return new_hashes

    def _process_query(self, conditions: Dict[str, str], depth=0) -> List[str]:

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
        raise NotImplementedError

    def query(
        self,
        query: str | dict,
        return_records=False
    ) -> GtfGff | List[dict]:

        if isinstance(query, str):
            query = self._process_query_string(query)

        hashes = self._process_query(query)
        records = self._get_records(hashes)

        if return_records:
            return records
        else:
            return self.gtf_gff_from_records(records)

    def export_records(self) -> List[dict]:
        # list comprehension ensures a copy of each record is made
        return [dict(record) for record in self.records.values()]

    def _remove_record(self, record_hash: str) -> None:

        record = self.records[record_hash]
        self.feature_index[record["feature"]].remove(record_hash)
        self.seqname_index[record["seqname"]].remove(record_hash)
        for key, value in record["attributes"].items():
            self.attribute_index[key][value].remove(record_hash)
        self._record_hashes.remove(record_hash)
        del self.records[record_hash]

    def remove_empty_field(self, field: str | tuple | list) -> None:
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

    # Define regex pattern for attributes
    attr_pattern = re.compile(r'(\w+)\s+"([^"]+)"')

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
            attributes = dict(attr_pattern.findall(fields[8]))

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
    for r in records:
        if r["start"] == r["end"]:
            continue
        interval_tree.add(
            Interval(
                r["start"],
                r["end"],
                data={k: v for k, v in r.items() if k not in ["start", "end"]}
            )
        )

    return interval_tree
