# file: scripts/generate_qc_report.py
import os

def get_html_list_of_files(species: str, files: list):
    html_list = ""
    for file in files:
        #get relative path
        file = os.path.relpath(file, os.path.join(species, "results", "qualitycontrol"))
        html_list += f"<li><a href='{file}'>{file}</a></li>"
    return html_list


def generate_qc_report(species: str, fastqc_raw, fastqc_trimmed, fastqc_quality_filtered,
                       multiqc_raw, multiqc_trimmed, multiqc_quality_filtered,
                       output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print(f"Generating quality check report for {species}")
    print(f"Output file: {output_file}")
    print(f"FastQC raw: {fastqc_raw}")
    print(f"FastQC trimmed: {fastqc_trimmed}")
    print(f"FastQC quality filtered: {fastqc_quality_filtered}")
    print(f"MultiQC raw: {multiqc_raw}")
    print(f"MultiQC trimmed: {multiqc_trimmed}")
    print(f"MultiQC quality filtered: {multiqc_quality_filtered}")

    with open(output_file, "w") as report:
        report.write(f"""
        <html>
            <head>
                <title>Quality Check Report: {species}</title>
            </head>
            <body>
                <h1>Quality Check Report: {species}</h1>

                <h2>MultiQC</h2>
                <h3>Raw Reads</h3>
                <ul>{get_html_list_of_files(species, [multiqc_raw])}</ul>
                <h3>Trimmed Reads</h3>
                <ul>{get_html_list_of_files(species, [multiqc_trimmed])}</ul>
                <h3>Quality Filtered Reads</h3>
                <ul>{get_html_list_of_files(species, [multiqc_quality_filtered])}</ul>

                <h2>FastQC</h2>
                <h3>Raw Reads</h3>
                <ul>{get_html_list_of_files(species, fastqc_raw)}</ul>
                <h3>Trimmed Reads</h3>
                <ul>{get_html_list_of_files(species, fastqc_trimmed)}</ul>
                <h3>Quality Filtered Reads</h3>
                <ul>{get_html_list_of_files(species, fastqc_quality_filtered)}</ul>

            </body>
        </html>
        """)


# --- Snakemake entrypoint ---
if __name__ == "__main__":
    species = snakemake.params.species

    generate_qc_report(
        species,
        fastqc_raw = list(snakemake.input.fastqc_raw),
        fastqc_trimmed = list(snakemake.input.fastqc_trimmed),
        fastqc_quality_filtered = list(snakemake.input.fastqc_quality_filtered),
        multiqc_raw = snakemake.input.multiqc_raw,
        multiqc_trimmed = snakemake.input.multiqc_trimmed,
        multiqc_quality_filtered = snakemake.input.multiqc_quality_filtered,
        output_file = str(snakemake.output[0])
    )
