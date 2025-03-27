import ftplib

# FTP server details
ftp_host = "ftp-private.ncbi.nlm.nih.gov"
ftp_user = "subftp"
ftp_pass = ""

# Your personal account folder as given by NCBI
account_folder = ""  # Replace this with your actual folder

# New subfolder for this submission (must be inside account folder)
submission_folder = "gorilla_upload_2025"

# Files to upload
files_to_upload = [
    "GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg-trimmed200.fa",
    "GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap2.p_ctg-trimmed200.fa",
    "GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.p_ctg-trimmed200.fa",
]

try:
    # Connect and log in
    ftp = ftplib.FTP(ftp_host)
    ftp.login(user=ftp_user, passwd=ftp_pass)
    print(f"Connected to {ftp_host}")

    # Go to your user account folder
    ftp.cwd(account_folder)
    print(f"Changed to account directory: {account_folder}")

    # Try to make the subfolder for this submission
    try:
        ftp.mkd(submission_folder)
        print(f"Created submission folder: {submission_folder}")
    except ftplib.error_perm as e:
        if "550" in str(e):
            print(f"Folder '{submission_folder}' may already exist.")
        else:
            raise

    # Now change into the new submission folder
    ftp.cwd(submission_folder)
    print(f"Changed into submission folder: {submission_folder}")

    # Upload the files
    for file_name in files_to_upload:
        with open(file_name, "rb") as file:
            print(f"Uploading {file_name}...")
            ftp.storbinary(f"STOR {file_name}", file)
            print(f"Uploaded: {file_name}")

    ftp.quit()
    print("All files uploaded and FTP connection closed.")

except ftplib.all_errors as e:
    print(f"FTP error: {e}")
