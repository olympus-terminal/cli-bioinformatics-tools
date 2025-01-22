import ftplib

# FTP server details
ftp_host = "ftp-private.ncbi.nlm.nih.gov"
ftp_user = "subftp"  # Replace with your username
ftp_pass = "your_password"  # Replace with your password

# List of files to upload
files_to_upload = [
    "file1.fq.tar.gz",
    "file2.fq.tar.gz",
    "file3.fq.tar.gz",
]

# Establish FTP connection
try:
    ftp = ftplib.FTP(ftp_host)
    ftp.login(user=ftp_user, passwd=ftp_pass)
    print(f"Connected to {ftp_host}")

    # Upload each file
    for file_name in files_to_upload:
        with open(file_name, "rb") as file:
            print(f"Uploading {file_name}...")
            ftp.storbinary(f"STOR /dir-to-upload-to/{file_name}", file)
            print(f"Successfully uploaded {file_name}")

    # Close the connection
    ftp.quit()
    print("All files uploaded and connection closed.")

except ftplib.all_errors as e:
    print(f"FTP error: {e}")
