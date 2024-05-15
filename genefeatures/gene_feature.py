from intervaltree import Interval, IntervalTree
from typing import Type, TypeVar
from .gtf_tools import GtfGff, parse_gtf
from .fasta_tools import extract_sequence

gtf = TypeVar("gtf", bound = "GtfGff")


class GeneFeature:


    def __init__(self, records: dict | list[dict] | gtf = None):

        self.locations = IntervalTree()
        self.transcript_ids = []
        self.gene_name = ""
        self.gene_id = ""

        if records is None:
            pass 
        elif isinstance(records, dict):
            self.add_record(records)
        else:
            self.add_records(records)

    
    def add_records(self, records: list[dict] | gtf):

        if isinstance(records, GtfGff):
            records = records.export_records()
        elif isinstance(records, list):
            pass
        else:
            raise TypeError(f"`records` argument must be a list[dict], or GtfGff; got type {type(records)}")

        for r in records:
            self.add_record(r)
    

    def add_record(self, r: dict):

        # ignore records without start and end
        try:
            loc = Interval(r.pop("start"), r.pop("end"), data = r)
            self.locations.add(loc)
        except KeyError:
            return

        if r["feature"] == "gene":
            if self.gene_id == "":
                self.gene_name = r["attributes"]["gene_name"]
                self.gene_id = r["attributes"]["gene_id"]
            else:
                # don't allow multiple genes in one GeneFeature object
                raise ValueError(
                    f"""
                    GeneFeature class only reperesents a single gene;
                    given {self.gene_name} and {r["attributes"]["gene_name"]}
                    """
                )

        elif r["feature"] == "transcript":
            tid = r["attributes"]["transcript_id"]
            if tid is not None and tid not in self.transcript_ids:
                self.transcript_ids.append(tid)


    def sort_locations(self):
        self.locations = IntervalTree(sorted(self.locations))
    

    @staticmethod
    def _get_interval_attr(interval: Interval, attribute: str | tuple | list):

        if isinstance(attribute, str): 
            return interval.data[attribute]

        elif isinstance(attribute, tuple | list):
            return interval.data[attribute[0]][attribute[1]]


    def partition_transcripts(self): 

        if len(self.transcript_ids) == 0:
            print("No transcripts found")
            return

        transcripts = {}
        for transcript_id in self.transcript_ids:
            transcripts[transcript_id] = IntervalTree() 

        for inter in self.locations:
            # skip if no transcript_id
            try:
                inter_tid = self._get_interval_attr(inter, ("attributes", "transcript_id"))
            except KeyError:
                continue

            # don't add whole gene or whole transcript intervals
            feature = self._get_interval_attr(inter, "feature")
            if feature in ("transcript", "gene"):
                continue

            # make new transcript entry if one doesn't exist
            if inter_tid not in transcripts.keys():
                self.transcript_ids.append(inter_tid)
                transcripts[inter_tid] = IntervalTree(inter)

            transcripts[inter_tid].add(inter)
        
        self.transcripts = transcripts


    def get_sequence(self, transcript_id):
        pass
    
    def translate_transcript(self, transcript_id):
        pass

    def mutate(self, transcript, level, change):
        pass
    
    @classmethod
    def from_fusion(cls):
        pass

