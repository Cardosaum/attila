#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pathlib
import subprocess
import time

files_path_raw = pathlib.Path("./data/raw/").glob("**/*")
files_path_parsed = pathlib.Path("./data/parsed/")
cdr_parser_path = pathlib.Path("./cdr.py").resolve()

for i in files_path_raw:
    if i.is_file() and str(i).endswith("txt") and ("aafreq" in str(i)):
        if "/VH/" in str(i):
            # we exclude the 2 firts elements because the directory
            # structure is data/raw/...
            # and we only need the data after /raw/ directory.
            print("=" * 30)
            time_start = time.perf_counter()
            filename_raw = "_".join(i.parts[3:])
            filename_parsed = pathlib.Path(pathlib.Path(filename_raw).stem + ".csv")
            file_path_parsed = pathlib.Path("./data/parsed/", filename_parsed).resolve()
            subprocess.run(["/usr/bin/python", str(cdr_parser_path), i.resolve(), file_path_parsed], check=True)
            time_stop = time.perf_counter() - time_start
            print(filename_raw)
            print(filename_parsed)
            print(f"Elapsed time: {time_stop}")

            print("=" * 30)
