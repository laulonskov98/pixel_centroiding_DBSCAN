import subprocess
import os

EXECUTABLE_PATH = "./dbscan_beta"

def compile_dbscan():
    # Specify the C++ source code file
    cpp_file = "./dbscan_beta.cpp"

    # Specify the output executable file name
    output_file = EXECUTABLE_PATH

    # Compile the C++ program using the C++ compiler (g++ in this example)
    try:
        compile_command = ["g++", cpp_file, "-o", output_file, "-O2"]
        subprocess.run(compile_command, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError as e:
        print("Compilation failed:", e)



def run_centroiding(input_path: str, output_path: str, correctiondata_path:str = None):
    """
    Description of run_centroiding.

    Args:
        input_path              (str):   "path/to/input.csv"
        output_path             (str):   "path/to/output.csv"
        arcorrectiondata_pathg3 (str):   "path/to/correctiondata.csv"

    Returns:
        None
    """

    file_path = EXECUTABLE_PATH

    if not os.path.exists(file_path):
        compile_dbscan()

    # Run the C++ program with arguments
    if correctiondata_path is not None:
        process = subprocess.Popen([EXECUTABLE_PATH, input_path, output_path, correctiondata_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    else:
        process = subprocess.Popen([EXECUTABLE_PATH, input_path, output_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    # Read the output and error (if any)
    output, error = process.communicate()
    if error != "":
        print(output, error)



