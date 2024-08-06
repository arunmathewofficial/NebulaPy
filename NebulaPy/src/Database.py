from urllib.parse import urljoin
import requests
import os
import re
import time
import sys
from requests.exceptions import Timeout
import version
import warnings
warnings.filterwarnings("ignore")


class download_database:

    def __init__(self, database_dir=None, verbose=False):
        self.database_dir = database_dir
        self.verbose = verbose
        #self.download_atlas_database()
        self.download_postdam_database()

    #########################################################################################
    # download progress bar
    #########################################################################################
    def download_progress(self, iteration, total_files, prefix='', fill='█', error=None):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total_files - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            fill       - Optional  : bar fill character (Str)
        """
        length = 20 # length of progress bar
        percent = ("{0:." + str(1) + "f}").format(100 * (iteration / float(total_files)))
        filled_length = int(length * iteration // total_files)
        bar = fill * filled_length + '-' * (length - filled_length)
        if error is None:
            sys.stdout.write(f'\r{prefix} downloading: |{bar}| {percent}% complete')
            sys.stdout.flush()
            # Print New Line on Complete
            if iteration == total_files:
                sys.stdout.write('\r' + ' ' * (len(f'{prefix} downloading: |{bar}| {percent}% complete')) + '\r')
                sys.stdout.flush()
                sys.stdout.write(f'\r{prefix} download completed\n')
        else:
            # Clear the progress bar line
            sys.stdout.write('\r' + ' ' * (len(f'{prefix} downloading: |{bar}| {percent}% complete')) + '\r')
            sys.stdout.flush()
            # Print the error message
            print(f'{error}')

    #########################################################################################
    # Verify atlas model in atlas database directory
    #########################################################################################
    def verify_atlas_model(self, model_dir, nfiles):
        """
        Verify whether the specified directory contains exactly `nfiles` and
        if all of those files are exactly 67 bytes in size.

        Args:
            model_dir (str): The path to the directory to be checked.
            nfiles (int): The expected number of files in the directory.

        Returns:
            bool: True if the directory contains exactly `nfiles` and all
             files are 69120 bytes in size, False otherwise.
        """
        try:
            files = os.listdir(model_dir)  # List all files in the specified directory
            if len(files) != nfiles:  # Check if the number of files matches the expected number
                return False  # Return False if the number of files does not match

            for file in files:  # Iterate over each file in the directory
                file_path = os.path.join(model_dir, file)  # Construct the full file path
                if os.path.isfile(file_path):  # Check if the current path is a file
                    if os.path.getsize(file_path) != 69120:  # Check if the file size is 69120 bytes
                        return False  # Return False if the file size is not 67 bytes
                else:
                    return False  # Return False if the current path is not a file

            return True  # Return True if all conditions are met
        except Exception as e:  # Catch any exceptions that occur
            print(f" error: {e}")  # Print the exception message
            return False  # Return False if an exception occurs

    #########################################################################################
    # Verify atlas model file in atlas database
    #########################################################################################
    def verify_atlas_model_file_exist(self, model_dir, file, filesize):
        """
        Check if a file exists in the specified directory and if
        its size match

        Parameters:
        model_dir (str): The directory to check.
        file (str): The file to look for.

        Returns:
        bool: True if the file exists and its size is 67 kB, False otherwise.
        """
        file_path = os.path.join(model_dir, file)
        if os.path.isfile(file_path):
            file_size = os.path.getsize(file_path)
            return file_size == filesize
        return False


    #########################################################################################
    # download atlas database
    #########################################################################################
    def download_atlas_database(self):
        atlas_base = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/ck04models/'
        atlas_models = ['ckm05/', 'ckm10/', 'ckm15/', 'ckm20/', 'ckm25/', 'ckp00/',
                        'ckp02/', 'ckp05/']
        Nfiles = [76, 76, 76, 76, 76, 76, 76, 76]
        atlas_info = ['AA_README', 'catalog.fits']
        atlas_dir = os.path.join(self.database_dir, 'Atlas/')
        os.makedirs(atlas_dir, exist_ok=True)

        # Downloading all atlas models ######################################################
        NModelExists = 0
        for model_index, model in enumerate(atlas_models):
            atlas_model_url = urljoin(atlas_base, model)
            model_name = model.replace('/', '')
            prefix_comment = f' Atlas model {model_name}'
            model_dir = os.path.join(atlas_dir, model)
            os.makedirs(model_dir, exist_ok=True)

            if self.verify_atlas_model(model_dir, Nfiles[model_index]):
                print(f" Atlas model {model_name} exists in database")
                NModelExists += 1
                continue

            else:
                response = requests.get(atlas_model_url)
                href_links = re.findall(r'href=[\'"]?([^\'" >]+)', response.text)
                # filter out non-FITS files
                href_links = [link for link in href_links if link.endswith('.fits')]

                # Download files with .fits extensions
                for file_index, href in enumerate(href_links, start=1):
                    # check if the file already exist in database, then skip
                    if self.verify_atlas_model_file_exist(model_dir, href, 69120):
                        continue
                    else:
                        file_url = os.path.join(atlas_model_url, href)
                        err = self.download_file(file_url, model_dir)
                        self.download_progress(file_index, Nfiles[model_index],
                                               prefix=prefix_comment, error=err)
                        if err is not None:
                            print(f"\033[91m error: PyMicroPION version"
                                  f" {version.__version__} database "
                                  f"download incomplete, please restart.\033[0m")
                            sys.exit(1)


        if NModelExists == 8:
            print(" all Atlas model exist in database")

        else:
            print(" Atlas model download completed successfully")
        # Downloading atlas info files ######################################################

    #########################################################################################
    # method to download single file
    #########################################################################################
    def download_file(self, url, dest_dir, timeout=10):
        try:
            # Create the destination directory if it does not exist
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)

            # Send a GET request to the URL with a specified timeout
            response = requests.get(url, stream=True, timeout=timeout)
            # Raise an error for bad HTTP status codes
            response.raise_for_status()

            # Extract the filename from the URL and create the full path
            filename = os.path.join(dest_dir, url.split('/')[-1])
            with open(filename, 'wb') as file:
                # Write the content to the file in chunks
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        file.write(chunk)
            return None

        except Timeout:
            # Handle the timeout exception and return an error message
            return " error: download timed out"
        except requests.RequestException as e:
            # Handle other HTTP exceptions and return an error message
            return f" error: {str(e)}"


    #########################################################################################
    # download postdam database
    #########################################################################################
    def download_postdam_database(self):
        import tarfile
        import time

        # Path to the tar.xz file
        tar_xz_path = 'hi.har.xz'
        # Directory where you want to extract the files
        extract_path = 'extracted_files'
        # Delay in seconds
        delay = 0.1  # Change this to the desired delay

        # Open the tar.xz file
        with tarfile.open(tar_xz_path, 'r:xz') as tar:
            # Iterate through the members of the tar archive
            for member in tar.getmembers():
                print(f'Extracting {member.name}')
                tar.extract(member, path=extract_path)
                time.sleep(delay)

        print(f'Files have been extracted to {extract_path}')



# Download all files from the directory
download_database('./PyMicroPION_database/')


