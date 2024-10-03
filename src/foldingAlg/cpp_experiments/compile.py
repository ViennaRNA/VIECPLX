#!/usr/bin/env python3

# 100% ChatGPT code...

import os
import subprocess

# Enable debug prints
DEBUG = True

def debug_print(message):
    if DEBUG:
        print(message)

# Function to retrieve default include paths from g++
def get_gpp_include_paths():
    try:
        # Run g++ to get the include paths it searches
        result = subprocess.run(
            ["g++", "-v", "-x", "c++", "-E", "-"], 
            input="", text=True, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL
        )
        output = result.stderr
        
        # Parse the output to find include paths
        include_paths = []
        start_collecting = False
        for line in output.splitlines():
            if line.startswith("#include <...> search starts here:"):
                start_collecting = True
                continue
            if start_collecting:
                if line.startswith("End of search list."):
                    break
                include_paths.append(line.strip())
        return include_paths
    except Exception as e:
        debug_print(f"Error while retrieving g++ include paths: {e}")
        return []

# Function to check if ViennaRNA's cofold.h exists in system include paths and return the include path if found
def find_viennarna_includes_in_system():
    # Check environment variables for additional include paths
    env_include_paths = []
    if "CPLUS_INCLUDE_PATH" in os.environ:
        env_include_paths.extend(os.environ["CPLUS_INCLUDE_PATH"].split(os.pathsep))
    if "CPATH" in os.environ:
        env_include_paths.extend(os.environ["CPATH"].split(os.pathsep))

    # Get default g++ include paths
    gpp_include_paths = get_gpp_include_paths()

    # Combine all include paths to search
    all_include_paths = gpp_include_paths + env_include_paths

    # Check if cofold.h exists in any of the paths
    for path in all_include_paths:
        cofold_path = os.path.join(path, "ViennaRNA", "cofold.h")
        if os.path.exists(cofold_path):
            debug_print(f"Found ViennaRNA header at: {cofold_path}")
            return path  # Return the include path where ViennaRNA was found

    return None  # Not found in system paths

# Function to check if ViennaRNA exists in a Conda environment and return the include path if found
def check_conda_viennarna():
    conda_env = os.getenv('CONDA_PREFIX')
    if conda_env:
        conda_include_path = os.path.join(conda_env, "include", "ViennaRNA", "cofold.h")
        if os.path.exists(conda_include_path):
            debug_print(f"Found ViennaRNA header in Conda environment at: {conda_include_path}")
            return os.path.join(conda_env, "include")  # Return the include path
    return None

# Function to compile the program
def compile_program(include_path=None, lib_path=None):
    gpp_command = ["g++", "beam_search_fold.cpp", "-o", "beam_search_fold", "-lRNA", "-lm", "-lstdc++", "-fopenmp", "-std=c++20"]
    
    # Add include and library paths if provided
    if include_path and lib_path:
        gpp_command.extend([f"-I{include_path}", f"-L{lib_path}"])

    # Print and run the compilation command
    print(f"Compiling with command: {' '.join(gpp_command)}")
    try:
        subprocess.run(gpp_command, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError:
        print("Error: Compilation failed.")

# Main function
def main():
    # Option 1: Check system include paths
    system_include_path = find_viennarna_includes_in_system()
    
    if system_include_path:
        print("Option 1: ViennaRNA found in system include paths.")
        compile_program(include_path=system_include_path, lib_path="/usr/lib")  # Assuming libraries are in /usr/lib
        return
    
    # Option 2: Check if Conda environment is active and ViennaRNA is found
    conda_include_path = check_conda_viennarna()
    
    if conda_include_path:
        print("Option 2: ViennaRNA found in Conda environment.")
        conda_lib_path = os.path.join(os.getenv('CONDA_PREFIX'), "lib")  # Conda libraries path
        compile_program(include_path=conda_include_path, lib_path=conda_lib_path)
        return
    
    # If neither option 1 nor option 2 is found
    print("Error: ViennaRNA header not found in system or Conda environment. Compilation aborted.")

# Run the main function
if __name__ == "__main__":
    main()