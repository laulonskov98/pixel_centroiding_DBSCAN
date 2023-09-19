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
        compile_command = ["g++", cpp_file, "-o", output_file, "-O3"]
        subprocess.run(compile_command, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError as e:
        print("Compilation failed:", e)



def run_centroiding(arg1, arg2, arg3=None):
    file_path = EXECUTABLE_PATH

    if not os.path.exists(file_path):
        compile_dbscan()

    # Run the C++ program with arguments
    if arg3 is not None:
        process = subprocess.Popen([EXECUTABLE_PATH, arg1, arg2, arg3], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    else:
        process = subprocess.Popen([EXECUTABLE_PATH, arg1, arg2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    # Read the output and error (if any)
    output, error = process.communicate()
    if error != "":
        print(output, error)



