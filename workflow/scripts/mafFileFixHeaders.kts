import java.io.*
import kotlin.system.exitProcess

// This script takes a multi alignment file (maf) created by the lastal alignment
// program and a second maf file created by the axtMaf program.  It grabs the header
// lines from the lastal maf, then adds them to the axtMaf file at the end of the
// existing headers in that file.
// Rational:  mafTools/mafCoverage needs the first line to be the version line as
// defined in the axtMaf file, but last-split needs the scoring matrix from the
// lastal maf file.
// For the gerp_align pipeline, last-split is not run until after chaining/netting,
// and is performed on the *.maf file created by the axtToMaf program.

// Check for parameters
if (!(args.contains("-lastalMaf")) || !(args.contains("-axtMaf")) || !(args.contains("-finalMaf"))) {
   println("must specify lastal maf with -lastalMaf parameter, axt maf with -axtMaf parameter  and output maf with -outputMaf parameter")
   exitProcess(1)
}

val lastalMaf = args[1 + args.indexOf("-lastalMaf")]
val axtMaf = args[1 + args.indexOf("-axtMaf")]
val finalMaf = args[1 + args.indexOf("-finalMaf")]

println("lastalMaf = $lastalMaf , axtMaf = $axtMaf , finalMaf = $finalMaf")

println("begin processing .... ")

// Step 1: get the header lines from the lastal maf file 

val lastalMafHeaders = getLastalHeaders(lastalMaf)

// Step2:  merge lastalMaf headers into axtToMaf file

mergeMafHeaders(axtMaf, lastalMafHeaders, finalMaf)

// This function parses the header lines from a maf (multiple alignment file)
// file and returns them in a List.
fun getLastalHeaders(mafFile: String ): List<String> {

    var headerList = mutableListOf<String>()
    val reader = BufferedReader(InputStreamReader(FileInputStream(mafFile)))
    var currentLine = reader.readLine()
    while(currentLine != null) {
        if(currentLine.startsWith("#")) {
            headerList.add(currentLine)
        } else {
            // done - don't read the rest of the file
            break
        }
        currentLine = reader.readLine()
    }
    reader.close()
    return headerList
}

// This method takes a list of MAF File header lines and inserts them
// after the existing header lines in the input Maf file.  A new file
// is created and written to "outputMaf".
// This new file has all the headers needed by different programs
// processing the MAF file in later steps of the pipeline.
// "last-split" needs headers from a maf file created by lastal
// "mafTools/mafCoverage" needs headers from a maf file created by axtToMaf
fun mergeMafHeaders(inputMaf: String, headerList: List<String>, outputMaf: String) {
    // This method takes the maf file created from axtToMaf and header lines from lastal's Maf
    // To the new file, add the header lines from the inputMaf file, then
    // add the lines from the headerList (header lines from lastal maf), then
    // add the non-header maf lines from the inputMaf file (axtToMaf file)

    val writer = BufferedWriter(OutputStreamWriter(FileOutputStream(outputMaf)))
    val reader = BufferedReader(InputStreamReader(FileInputStream(inputMaf)))

    // Grab and write the maf data lines (from the axtToMaf file)
    var currentLine = reader.readLine() // next header line
    var lastalHeadersAdded = false
    while(currentLine != null) {
        // First: write the axtToMaf header lines - because some pograms e.g. mafTools/mafCoverage
        // need the first line to be a version line in axttoMaf format
        // Second: write the lastal maf headers lines:  because last-split needs the scoring matrix
        // Third: write the non-header maf lines
        if(!(currentLine.startsWith("#")) && lastalHeadersAdded == false) {
            for (line in headerList) {
                writer.write(line)
                writer.write("\n")
            }
            lastalHeadersAdded = true
        } else {
            writer.write(currentLine)
            writer.write("\n")

        } 
        currentLine = reader.readLine()
    }
    reader.close()
    writer.close()
}
println("done! writing to $finalMaf !!")
