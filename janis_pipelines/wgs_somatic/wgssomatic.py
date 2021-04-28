from janis_bioinformatics.data_types import FastqGzPair
from janis_core import String, Array, InputDocumentation, InputQualityType

from janis_pipelines.wgs_somatic.wgssomatic_variantsonly import (
    WGSSomaticMultiCallersVariantsOnly,
)
from janis_pipelines.wgs_somatic_gatk.wgssomaticgatk import WGSSomaticGATK


class WGSSomaticMultiCallers(WGSSomaticMultiCallersVariantsOnly, WGSSomaticGATK):
    def id(self):
        return "WGSSomaticMultiCallers"

    def friendly_name(self):
        return "WGS Somatic (Multi callers)"

    def version(self):
        return "1.4.0"

    def constructor(self):
        # don't call super()

        self.add_inputs()
        self.add_preprocessing_steps()
        self.add_gridss(
            normal_bam_source=self.normal.out_bam, tumor_bam_source=self.tumor.out_bam
        )
        self.add_facets(
            normal_bam_source=self.normal.out_bam, tumor_bam_source=self.tumor.out_bam
        )
        self.add_gatk_variantcaller(
            normal_bam_source=self.normal.out_bam, tumor_bam_source=self.tumor.out_bam
        )
        self.add_vardict_variantcaller(
            normal_bam_source=self.normal.out_bam, tumor_bam_source=self.tumor.out_bam
        )
        self.add_strelka_variantcaller(
            normal_bam_source=self.normal.out_bam, tumor_bam_source=self.tumor.out_bam
        )
        self.add_combine_variants(
            normal_bam_source=self.normal.out_bam, tumor_bam_source=self.tumor.out_bam
        )

    def add_inputs(self):
        # INPUTS
        self.input(
            "normal_inputs",
            Array(FastqGzPair),
            doc=InputDocumentation(
                "An array of NORMAL FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads",
                quality=InputQualityType.user,
                example='["normal_R1.fastq.gz", "normal_R2.fastq.gz"]',
            ),
        )
        self.input(
            "tumor_inputs",
            Array(FastqGzPair),
            doc=InputDocumentation(
                "An array of TUMOR FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads",
                quality=InputQualityType.user,
                example='["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"]',
            ),
        )
        self.input(
            "normal_name",
            String(),
            doc=InputDocumentation(
                "Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem",
                quality=InputQualityType.user,
                example="NA24385_normal",
            ),
        )
        self.input(
            "tumor_name",
            String(),
            doc=InputDocumentation(
                "Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem",
                quality=InputQualityType.user,
                example="NA24385_tumor",
            ),
        )

        self.add_inputs_for_reference()
        self.add_inputs_for_intervals()
        self.add_inputs_for_configuration()


if __name__ == "__main__":
    import os.path

    w = WGSSomaticMultiCallers()
    args = {
        "to_console": False,
        "to_disk": True,
        "validate": True,
        "export_path": os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "{language}"
        ),
    }
    w.translate("cwl", **args)
    w.translate("wdl", **args)
