''' The following scipt takes in a directory path as input, and 
    unzips all the gzip and zip files in that directory.
    It will not overwrite existing files.
    Usage: python unzip_if_needed.py <directory_path>'''
import os
import sys
import gzip
import zipfile

def ungzip_file(gzip_path, dest_path):
    with gzip.open(gzip_path, 'rb') as f_in, open(dest_path, 'wb') as f_out:
        f_out.write(f_in.read())
    print(f"Extracted gzip: {gzip_path} -> {dest_path}")

def unzip_file(zip_path, dest_dir):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for member in zip_ref.namelist():
            target_path = os.path.join(dest_dir, member)
            if not os.path.exists(target_path):
                zip_ref.extract(member, dest_dir)
                print(f"Extracted zip: {member} -> {target_path}")
            else:
                print(f"Skipped existing file: {target_path}")

def main(path):
    if not os.path.isdir(path):
        print(f"Invalid path: {path}")
        return

    for root, _, files in os.walk(path):
        for filename in files:
            full_path = os.path.join(root, filename)

            if filename.endswith('.gzip'):
                dest_path = full_path[:-5]  # remove .gzip
                if not os.path.exists(dest_path):
                    ungzip_file(full_path, dest_path)
                else:
                    print(f"Skipped existing file: {dest_path}")

            elif filename.endswith('.zip'):
                unzip_file(full_path, root)

            if filename.endswith('.gz'):
                dest_path = full_path[:-3]  # remove .gz
                if not os.path.exists(dest_path):
                    ungzip_file(full_path, dest_path)
                else:
                    print(f"Skipped existing file: {dest_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python unzip_if_needed.py <directory_path>")
    else:
        main(sys.argv[1])