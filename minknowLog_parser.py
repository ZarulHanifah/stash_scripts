#!/usr/bin/env python

import re
import argparse
from typing import List, Dict, Optional


class Record:
    def __init__(self, date: str, time: str, info: str, attrs: Optional[Dict[str, str]] = None):
        self.date = date
        self.time = time
        self.info = info
        self.attrs = attrs or {}

    def __repr__(self):
        return f"<Record {self.date} {self.time} {self.info}>"


class LogParser:
    RECORD_REGEX = re.compile(r"^(\d{4}-\d{2}-\d{2})\s+(\d{2}:\d{2}:\d{2}\.\d+)\s+(INFO:.*)")

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.records: List[Record] = []

    def parse(self):
        with open(self.filepath, "r", encoding="utf-8") as f:
            current_record: Optional[Record] = None
            for line in f:
                line = line.rstrip("\n")
                match = self.RECORD_REGEX.match(line)
                if match:
                    if current_record:
                        self.records.append(current_record)
                    date, time, info = match.groups()
                    info = re.sub("INFO: \(\s*\) \(", "\1", info)
                    # info = re.sub("INFO: \(.*\) attrs.*", "\1", info)
                    info_match = re.search(r'INFO:\s+(\w+)', info)
                    info = info_match.group(1)
                    current_record = Record(date, time, info)
                elif current_record and line.startswith("    "):
                    key_value = line.strip().split(": ", 1)
                    if len(key_value) == 2:
                        key, value = key_value
                        current_record.attrs[key] = value
            if current_record:
                self.records.append(current_record)

    def get_records(self) -> List[Record]:
        return self.records


def main():
    parser = argparse.ArgumentParser(description="Parse a log file into structured records.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input log file")
    parser.add_argument("-t", "--target", default=None, help="Find target")
    args = parser.parse_args()

    log_parser = LogParser(args.input)
    log_parser.parse()

    for record in log_parser.get_records():
    # for record in log_parser.get_records()[:30]:
        if args.target != None:
            target = args.target
            if record.info in target:
                print(record)
                for k, v in record.attrs.items():
                    print(f"    {k}: {v}")
        else:
            print(record)
            for k, v in record.attrs.items():
                print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
