import os
import pandas
import ipaddress

'''
Script to process a Caida2019 network trace to csv in the desired format.
Requirements: gzip, tcpdump, tshark -> installed via $ sudo apt-get install gzip tcpdump tshark

Given a network trace, such as equinix-nyc.dirA.yyyymmdd-hhmmss.UTC.anon.pcap.gz
1. It is unzipped by gzip to a obtain a pcap
2. The pcap is split into smaller pcaps of 100MB each, as tshark cannot process the large pcaps
3. Each pcap split is processed by tshark which outputs a csv with [timestamp, ip.src, ip.dst]

'''

CAIDA_DIR = "traces"   # folder with caida traces
REMOVE_PROCESSED_FILES = True  # let the script remove processed files to help save space
SPLIT_SIZE = 200  # MB size of caida trace splits
FIX_FORMAT = True  # fix format of [timestamp (datetime), ip.src (str), ip.dst (str)] -> [timestamp (bigint), ip.src (int), ip.dst (int)]
                    # i.e., [Jan 17, 2019 14:10:00.000003000 CET, 46.244.245.152, 171.198.20.110] -> [1547728200000003000, 778246152, 2886798318]

CREATE_FOLDERS = False  # let scrip create split_folder and output_folder
split_folder = CAIDA_DIR + "_splits"
output_folder = CAIDA_DIR + "_output"


def fix_format(df):
    renaming = {'frame.time':'timestamp', 'ip.src':'ip.src', 'ip.dst':'ip.dst'}
    desired_types = {'timestamp':'string', 'ip.src':'string', 'ip.dst':'string'}
    
    df = df.rename(columns=renaming).astype(desired_types)
    df = df.dropna()
    df['timestamp'] = 0 #(pandas.to_datetime(df['timestamp']).astype(int) / 1000).astype(int)  # note, this step is slow, if timestamp is not required this can be set to any number
    df['ip.src'] = df['ip.src'].apply(lambda x: int(ipaddress.ip_address(x)))
    df['ip.dst'] = df['ip.dst'].apply(lambda x: int(ipaddress.ip_address(x)))


    return df

# Create relevant folders
if CREATE_FOLDERS:
    if (os.system("mkdir {} && mkdir {}".format(split_folder, output_folder)) != 0):  # folder to put splits in
        print("Error creating splits and output folder")
        exit()

# For all files in caida_dir
for file in sorted(os.listdir(CAIDA_DIR)):
    

    print("Unzipping file", file)
    if os.system("gzip -dc {} > {}".format(os.path.join(CAIDA_DIR, file), os.path.join(CAIDA_DIR, file.split(".gz")[0]))) != 0:
        print("Error unzipping")
        break

    print("Splitting file", file, "into", SPLIT_SIZE, "MB chunks")
    if os.system("tcpdump -r {} -w {} -C {}".format(os.path.join(CAIDA_DIR, file.split(".gz")[0]), os.path.join(split_folder, file.split(".pcap")[0]), SPLIT_SIZE)) != 0:
        print("Error splitting")
        break
    elif REMOVE_PROCESSED_FILES:
        os.system("rm {}".format(os.path.join(CAIDA_DIR, file.split(".gz")[0])))

    # Process splits
    split_success = True
    for split_file in sorted(os.listdir(split_folder)):
        print("\tProcessing split", split_file)
        if (os.system("tshark -r {} -T fields -e frame.time -e ip.src -e ip.dst -E header=y -E separator=, -E quote=d -E occurrence=f > {}/{}.csv".format(os.path.join(split_folder, split_file), output_folder, split_file))) != 0:
            print("Tshark process failed")
            split_success = False
            break

        if FIX_FORMAT:
            fix_format(pandas.read_csv(os.path.join(output_folder, split_file+".csv"))).to_csv(os.path.join(output_folder, split_file+".csv"), index=False)

        # if this statement is reached, tshark succeeded, thus remove spit if flag is set
        if REMOVE_PROCESSED_FILES:
            print("\tRemoving split")
            os.system("rm {}".format(os.path.join(split_folder, split_file)))

    if not split_success:
        break
    
    # Remove original file is flag set
    if REMOVE_PROCESSED_FILES:
        print("Removing trace", file)
        os.system("rm {}".format(os.path.join(CAIDA_DIR, file)))

if CREATE_FOLDERS:
    print("Attempting to remove split folder")
    os.system("rm -r {}".format(split_folder))