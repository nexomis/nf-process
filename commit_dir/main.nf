// Groovy function to determine the output directory name
// Defined globally in the script file to be accessible in output directive
def getOutputDirName(inputPath) {
    out_dir="${inputPath}"
    if (inputPath.isDirectory()) {
        // If it's a directory, the output name *is* the directory name
        return out_dir
    } else  {
        // If it's a file, derive name by removing extensions
        out_dir = out_dir.replaceAll('(?i)\\.gz$', '')
        out_dir = out_dir.replaceAll('(?i)\\.bz2$', '')
        out_dir = out_dir.replaceAll('(?i)\\.z$', '')
        out_dir = out_dir.replaceAll('(?i)\\.bz$', '')
        out_dir = out_dir.replaceAll('(?i)\\.tar$', '')

        // Check if name changed significantly or if it was not a recognized archive
        if (out_dir=="${inputPath}") {
             out_dir = "${out_dir}_extracted"
        }
        // Handle cases where extension removal leaves an empty string
        if (out_dir.isEmpty()) {
            out_dir = "extracted"
        }
        return out_dir
    }
}

process COMMIT_DIR {
    tag "$meta.id"
    label 'process_low'
    container "python:3.14.0a7-bookworm" // Or a suitable image with file, tar

    input:
    tuple val(meta), path(input_item)

    output:
    // Use the globally defined function here
    tuple val(meta), path("${getOutputDirName(input_item)}"), emit: committed_dir

    script:
    // Call the function again inside the script context to get the name for bash
    def out_dir = getOutputDirName(input_item)
    """
    #!/bin/bash
    set -e
    echo "Processing input: ${input_item}"

    if [ -d "${input_item}" ]; then
        echo "Input is a directory. Output directive will handle it."
    elif [ -f "${input_item}" ]; then
        echo "Input is a file. Checking type and extracting to '${out_dir}'."
        mime_type=\$(file -L -b --mime-type "${input_item}")
        echo "Detected MIME type: \$mime_type"
        mkdir "${out_dir}"

        case "\$mime_type" in
            application/gzip|application/x-gzip)
                echo "Extracting Gzip archive..."
                tar -xzf "${input_item}" -C "${out_dir}" || { echo "Error extracting Gzip archive"; exit 1; }
                ;;
            application/x-bzip2)
                echo "Extracting Bzip2 archive..."
                tar -xjf "${input_item}" -C "${out_dir}" || { echo "Error extracting Bzip2 archive"; exit 1; }
                ;;
            application/x-tar)
                echo "Extracting Tar archive..."
                tar -xf "${input_item}" -C "${out_dir}" || { echo "Error extracting Tar archive"; exit 1; }
                ;;
            *)
                echo "Error: Unsupported file type '\$mime_type' for input ${input_item}"
                exit 1
                ;;
        esac
        echo "Extraction complete."
    else
        echo "Error: Input ${input_item} is neither a directory nor a regular file."
        exit 1
    fi

    echo "Output directory prepared."
    """

    stub:
    // Stub needs to reflect the output structure, using the global function
    def out_dir_stub = getOutputDirName(input_item)
    """
    mkdir ${out_dir_stub}
    touch ${out_dir_stub}/stub_file
    """
}
