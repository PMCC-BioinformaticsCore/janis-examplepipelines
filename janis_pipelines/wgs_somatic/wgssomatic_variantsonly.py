from datetime import date

from janis_bioinformatics.data_types import Bed, BedTabix, Vcf, CompressedVcf
from janis_bioinformatics.tools.bcftools import BcfToolsSort_1_9
from janis_bioinformatics.tools.common import GATKBaseRecalBQSRWorkflow_4_1_3
from janis_bioinformatics.tools.gatk4 import Gatk4GatherVcfs_4_1_3
from janis_bioinformatics.tools.htslib import BGZipLatest
from janis_bioinformatics.tools.papenfuss import Gridss_2_6_2
from janis_bioinformatics.tools.pmac import (
    CombineVariants_0_0_8,
    GenerateVardictHeaderLines,
    AddBamStatsSomatic_0_1_0,
    GenerateIntervalsByChromosome,
)
from janis_bioinformatics.tools.variantcallers import GatkSomaticVariantCaller_4_1_3
from janis_bioinformatics.tools.variantcallers.illuminasomatic_strelka import (
    IlluminaSomaticVariantCaller,
)
from janis_bioinformatics.tools.variantcallers.vardictsomatic_variants import (
    VardictSomaticVariantCaller,
)
from janis_bioinformatics.tools.common.facetsWorkflow import FacestWorkflow
from janis_core import (
    Array,
    Float,
    Int,
    Boolean,
    String,
    WorkflowMetadata,
    InputDocumentation,
    InputQualityType,
    StringFormatter,
)
from janis_core.operators.standard import FirstOperator
from janis_unix.tools import UncompressArchive

from janis_pipelines.wgs_somatic_gatk.wgssomaticgatk_variantsonly import (
    WGSSomaticGATKVariantsOnly,
    INPUT_DOCS,
)


class WGSSomaticMultiCallersVariantsOnly(WGSSomaticGATKVariantsOnly):
    def id(self):
        return "WGSSomaticMultiCallersVariantsOnly"

    def friendly_name(self):
        return "WGS Somatic (Multi callers) [VARIANTS only]"

    def version(self):
        return "1.4.0"

    def constructor(self):
        # don't call super()

        self.add_inputs()
        self.add_gridss(
            normal_bam_source=self.normal_bam, tumor_bam_source=self.tumor_bam
        )
        self.add_facets(
            normal_bam_source=self.normal_bam, tumor_bam_source=self.tumor_bam
        )
        self.add_gatk_variantcaller(
            normal_bam_source=self.normal_bam, tumor_bam_source=self.tumor_bam
        )
        self.add_vardict_variantcaller(
            normal_bam_source=self.normal_bam, tumor_bam_source=self.tumor_bam
        )
        self.add_strelka_variantcaller(
            normal_bam_source=self.normal_bam, tumor_bam_source=self.tumor_bam
        )
        self.add_combine_variants(
            normal_bam_source=self.normal_bam, tumor_bam_source=self.tumor_bam
        )

    def add_inputs_for_intervals(self):
        super().add_inputs_for_intervals()

        self.input("vardict_intervals", Array(Bed), doc=INPUT_DOCS["vardict_intervals"])
        self.input("strelka_intervals", BedTabix, doc=INPUT_DOCS["strelka_intervals"])

    def add_inputs_for_configuration(self):
        super().add_inputs_for_configuration()
        # facets
        self.input("pseudo_snps", Int(optional=True))
        self.input("max_depth", Int(optional=True))
        self.input("everything", Boolean(optional=True))
        self.input("genome", String(optional=True))
        self.input("cval", Int(optional=True))
        self.input("purity_cval", Int(optional=True))
        self.input("normal_depth", Int(optional=True))
        # vardict
        self.input(
            "allele_freq_threshold",
            Float,
            default=0.05,
            doc=InputDocumentation(
                "The threshold for VarDict's allele frequency, default: 0.05 or 5%",
                quality=InputQualityType.configuration,
            ),
        )
        self.input("minMappingQual", Int(optional=True))
        self.input("filter", String(optional=True))

    def add_gridss(self, normal_bam_source, tumor_bam_source):

        # GRIDSS
        self.step(
            "vc_gridss",
            Gridss_2_6_2(
                bams=[normal_bam_source, tumor_bam_source],
                reference=self.reference,
                blacklist=self.gridss_blacklist,
            ),
        )

        # GRIDSS
        self.output(
            "out_gridss_assembly",
            source=self.vc_gridss.assembly,
            output_folder=[
                "gridss",
            ],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_gridss",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Assembly returned by GRIDSS",
        )
        self.output(
            "out_variants_gridss",
            source=self.vc_gridss.out,
            output_folder=[
                "gridss",
            ],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_gridss",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Variants from the GRIDSS variant caller",
        )

    def add_facets(self, normal_bam_source, tumor_bam_source):

        self.step(
            "vc_facets",
            FacestWorkflow(
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
                snps_dbsnp=self.snps_dbsnp,
                pseudo_snps=self.pseudo_snps,
                max_depth=self.max_depth,
                everything=self.everything,
                genome=self.genome,
                cval=self.cval,
                purity_cval=self.purity_cval,
                normal_depth=self.normal_depth,
            ),
        )
        self.output(
            "out_facets_summary",
            source=self.vc_facets.out_summary,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}.txt",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        )
        self.output(
            "out_facets_purity_png",
            source=self.vc_facets.out_purity_png,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}_purity.png",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_purity_seg",
            source=self.vc_facets.out_purity_seg,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}_purity.seg",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_purity_rds",
            source=self.vc_facets.out_purity_rds,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}_purity.rds",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_hisens_png",
            source=self.vc_facets.out_hisens_png,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}_hisens.png",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_hisens_seg",
            source=self.vc_facets.out_hisens_seg,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}_hisens.seg",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_hisens_rds",
            source=self.vc_facets.out_hisens_rds,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}_hisens.rds",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_arm_level",
            source=self.vc_facets.out_arm_level,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}.arm_level.txt",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_gene_level",
            source=self.vc_facets.out_gene_level,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}.gene_level.txt",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),
        self.output(
            "out_facets_qc",
            source=self.vc_facets.out_qc,
            output_folder="facets",
            output_name=StringFormatter(
                "{tumour}--{normal}.qc.txt",
                tumour=self.tumor_name,
                normal=self.normal_name,
            ),
        ),

    def add_gatk_variantcaller(self, normal_bam_source, tumor_bam_source):
        """
        Reimplemented because need steps for combine
        """

        if "generate_gatk_intervals" in self.step_nodes:
            generated_intervals = self.generate_gatk_intervals.out_regions
        else:
            generated_intervals = self.step(
                "generate_gatk_intervals",
                GenerateIntervalsByChromosome(reference=self.reference),
                when=self.gatk_intervals.is_null(),
            ).out_regions

        intervals = FirstOperator([self.gatk_intervals, generated_intervals])

        recal_ins = {
            "reference": self.reference,
            "intervals": intervals,
            "snps_dbsnp": self.snps_dbsnp,
            "snps_1000gp": self.snps_1000gp,
            "known_indels": self.known_indels,
            "mills_indels": self.mills_indels,
        }

        self.step(
            "bqsr_normal",
            GATKBaseRecalBQSRWorkflow_4_1_3(bam=normal_bam_source, **recal_ins),
            scatter="intervals",
        )

        self.step(
            "bqsr_tumor",
            GATKBaseRecalBQSRWorkflow_4_1_3(bam=tumor_bam_source, **recal_ins),
            scatter="intervals",
        )

        self.step(
            "vc_gatk",
            GatkSomaticVariantCaller_4_1_3(
                normal_bam=self.bqsr_normal.out,
                tumor_bam=self.bqsr_tumor.out,
                normal_name=self.normal_name,
                intervals=intervals,
                reference=self.reference,
                gnomad=self.gnomad,
                panel_of_normals=self.panel_of_normals,
            ),
            scatter=["intervals", "normal_bam", "tumor_bam"],
        )

        self.step("vc_gatk_merge", Gatk4GatherVcfs_4_1_3(vcfs=self.vc_gatk.out))
        self.step("vc_gatk_compress_for_sort", BGZipLatest(file=self.vc_gatk_merge.out))
        self.step(
            "vc_gatk_sort_combined",
            BcfToolsSort_1_9(
                vcf=self.vc_gatk_compress_for_sort.out.as_type(CompressedVcf)
            ),
        )
        self.step(
            "vc_gatk_uncompress_for_combine",
            UncompressArchive(file=self.vc_gatk_sort_combined.out),
        )

        # VCF
        self.output(
            "out_variants_gatk",
            source=self.vc_gatk_sort_combined.out,
            output_folder=[
                "vcf",
            ],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_gatk",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Merged variants from the GATK caller",
        )
        self.output(
            "out_variants_split",
            source=self.vc_gatk.out,
            output_folder=[
                "vcf",
                "GATKByInterval",
            ],
            doc="Unmerged variants from the GATK caller (by interval)",
        )

    def add_strelka_variantcaller(self, normal_bam_source, tumor_bam_source):
        self.step(
            "vc_strelka",
            IlluminaSomaticVariantCaller(
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                intervals=self.strelka_intervals,
                reference=self.reference,
            ),
        )
        self.step("vc_strelka_compress", BGZipLatest(file=self.vc_strelka.out))

        self.output(
            "out_variants_strelka",
            source=self.vc_strelka_compress.out.as_type(CompressedVcf),
            output_folder=[
                "vcf",
            ],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_strelka",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Variants from the Strelka variant caller",
        )

    def add_vardict_variantcaller(self, normal_bam_source, tumor_bam_source):
        self.step(
            "generate_vardict_headerlines",
            GenerateVardictHeaderLines(reference=self.reference),
        )
        self.step(
            "vc_vardict",
            VardictSomaticVariantCaller(
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                normal_name=self.normal_name,
                tumor_name=self.tumor_name,
                header_lines=self.generate_vardict_headerlines.out,
                intervals=self.vardict_intervals,
                reference=self.reference,
                allele_freq_threshold=self.allele_freq_threshold,
                minMappingQual=self.minMappingQual,
                filter=self.filter,
            ),
            scatter="intervals",
        )
        self.step("vc_vardict_merge", Gatk4GatherVcfs_4_1_3(vcfs=self.vc_vardict.out))
        self.step(
            "vc_vardict_compress_for_sort", BGZipLatest(file=self.vc_vardict_merge.out)
        )
        self.step(
            "vc_vardict_sort_combined",
            BcfToolsSort_1_9(
                vcf=self.vc_vardict_compress_for_sort.out.as_type(CompressedVcf)
            ),
        )
        self.step(
            "vc_vardict_uncompress_for_combine",
            UncompressArchive(file=self.vc_vardict_sort_combined.out),
        )

        self.output(
            "out_variants_vardict",
            source=self.vc_vardict_sort_combined.out,
            output_folder=[
                "vcf",
            ],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_vardict",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Merged variants from the VarDict caller",
        )
        self.output(
            "out_variants_vardict_split",
            source=self.vc_vardict.out,
            output_folder=[
                "vcf",
                "VardictByInterval",
            ],
            doc="Unmerged variants from the GATK caller (by interval)",
        )

    def add_combine_variants(self, normal_bam_source, tumor_bam_source):
        self.step(
            "combine_variants",
            CombineVariants_0_0_8(
                normal=self.normal_name,
                tumor=self.tumor_name,
                vcfs=[
                    self.vc_gatk_uncompress_for_combine.out.as_type(Vcf),
                    self.vc_strelka.out,
                    self.vc_vardict_uncompress_for_combine.out.as_type(Vcf),
                ],
                type="somatic",
                columns=["AD", "DP", "GT"],
            ),
        )

        self.step("combined_compress", BGZipLatest(file=self.combine_variants.out))
        self.step(
            "combined_sort",
            BcfToolsSort_1_9(vcf=self.combined_compress.out.as_type(CompressedVcf)),
        )
        self.step("combined_uncompress", UncompressArchive(file=self.combined_sort.out))

        self.step(
            "combined_addbamstats",
            AddBamStatsSomatic_0_1_0(
                normal_id=self.normal_name,
                tumor_id=self.tumor_name,
                normal_bam=normal_bam_source,
                tumor_bam=tumor_bam_source,
                vcf=self.combined_uncompress.out.as_type(Vcf),
                reference=self.reference,
            ),
        )

        self.output(
            "out_variants",
            source=self.combined_addbamstats.out,
            output_folder=[
                "vcf",
            ],
            output_name=StringFormatter(
                "{tumor_name}--{normal_name}_combined",
                tumor_name=self.tumor_name,
                normal_name=self.normal_name,
            ),
            doc="Combined variants from GATK, VarDict and Strelka callers",
        )

    def bind_metadata(self):
        meta: WorkflowMetadata = super().bind_metadata() or self.metadata

        meta.keywords = [
            "wgs",
            "cancer",
            "somatic",
            "variants",
            "gatk",
            "vardict",
            "strelka",
            "gridss",
        ]
        meta.contributors = ["Michael Franklin", "Richard Lupat", "Jiaan Yu"]
        meta.dateCreated = date(2018, 12, 24)
        meta.dateUpdated = date(2020, 8, 19)
        meta.short_documentation = "A somatic tumor-normal variant-calling WGS pipeline using GATK, VarDict and Strelka2."
        meta.documentation = """\
This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs:

- Takes raw sequence data in the FASTQ format;
- align to the reference genome using BWA MEM;
- Marks duplicates using Picard;
- Call the appropriate somatic variant callers (GATK / Strelka / VarDict);
- Outputs the final variants in the VCF format.

**Resources**

This pipeline has been tested using the HG38 reference set, available on Google Cloud Storage through:

- https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

This pipeline expects the assembly references to be as they appear in that storage \
    (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
The known sites (snps_dbsnp, snps_1000gp, known_indels, mills_indels) should be gzipped and tabix indexed.
"""


if __name__ == "__main__":
    import os.path

    w = WGSSomaticMultiCallersVariantsOnly()
    args = {
        "to_console": False,
        "to_disk": False,
        "validate": True,
        "export_path": os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "{language}"
        ),
    }
    # w.translate("cwl", **args)
    w.translate("wdl", **args)
    # WGSSomaticMultiCallersVariantsOnly().translate("wdl")
