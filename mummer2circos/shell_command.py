import subprocess

def shell_command(command):
    """execute shell command and return std output and std error"""
    process = subprocess.Popen(command, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True)
    stdout_str, stderr_str = process.communicate()
    return(stdout_str, stderr_str, process.returncode)
