// Helper function to parse an argument string into a map of options and a list of flags
def parseArguments(String argString) {
    def options = [:]
    def flags = []
    if (argString == null || argString.trim().isEmpty()) {
        return [options: options, flags: flags]
    }
    // Regex to split by space, but respect quotes. Handles simple cases.
    // Allows for quoted arguments, e.g. -a "adapter seq"
    def tokens = argString.tokenize(/(?<!\\)".*?(?<!\\)"|\S+/)
                        .collect { it -> it.startsWith('"') && it.endsWith('"') ? it.substring(1, it.length() - 1) : it }

    String pendingFlag = ""
    tokens.each { token ->
        if (token.startsWith("-")) {
            // Current token is a flag - handle any pending flag first
            if (!pendingFlag.isEmpty()) {
                // Previous flag had no value, so it's a standalone flag
                flags.add(pendingFlag)
            }
            pendingFlag = token
        } else {
            if (!pendingFlag.isEmpty()) {
                // This token is the value for the pending flag
                options[pendingFlag] = token
                pendingFlag = ""
            }
            // If no pending flag, this token is ignored (shouldn't happen in well-formed args)
        }
    }
    
    if (!pendingFlag.isEmpty()) {
        flags = flags + [pendingFlag]
    }
    return [options: options, flags: flags]
}

// Function to combine default and meta arguments
// Meta arguments (from samplesheet) take precedence over default arguments (from ext.config)
def combineArgs(String defaultArgStr, String metaArgStr) {
    def parsedDefaults = parseArguments(defaultArgStr ?: "")
    def parsedMeta = parseArguments(metaArgStr ?: "")

    // Merge options, meta takes precedence
    def finalOptions = parsedDefaults.options
    finalOptions.putAll(parsedMeta.options) // Meta overrides defaults for same option keys

    // Merge flags, ensuring uniqueness
    def finalFlags = parsedDefaults.flags
    finalFlags.addAll(parsedMeta.flags)

    // Reconstruct the argument string
    def result = []
    finalOptions.each { key, value -> result.add(key); result.add("\"${value}\"") } // Quote values
    finalFlags.each { flag -> result.add(flag) }
    
    return result.join(' ').trim()
}

process CUTADAPT {
    tag "${meta.id}"
    cpus 4
    memory 8.GB

    container "quay.io/biocontainers/cutadapt:5.0--py312h0fa9677_0"

    input:
    tuple val(meta), path(reads, arity: 1..2, stageAs: 'input_raw/*')

    output:
    tuple val(meta), path("*.fq.gz", arity: '1..2'), emit: reads
    tuple val(meta), path("*.cutadapt.json"), emit: report_json

    script:
    def is_paired = reads.size() > 1
    def prefix = meta.id // Use meta.id for output prefix as per convention

    def combined_args_str = combineArgs(task.ext.args ?: "", meta.args_cutadapt ?: "")

    def input_files_str = is_paired ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    def output_files_option_str = is_paired ? "-o ${prefix}_R1.fq.gz -p ${prefix}_R2.fq.gz" : "-o ${prefix}_R1.fq.gz"
    
    """
    cutadapt \\
        -j $task.cpus \\
        $combined_args_str \\
        --json ${prefix}.cutadapt.json \\
        $output_files_option_str \\
        $input_files_str

    """

    stub:
    def prefix = meta.id
    def is_paired = reads.size() > 1
    def touch_files = is_paired ? "${prefix}_R1.fq.gz ${prefix}_R2.fq.gz" : "${prefix}_R1.fq.gz"
    """
    touch $touch_files
    touch ${prefix}.cutadapt.json
    """
}
