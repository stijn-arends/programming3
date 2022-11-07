"""
module
"""

import csv
from pathlib import Path

import pyspark.sql.functions as f
from pyspark import SparkContext
from pyspark.sql import SparkSession

# from pyspark.sql.functions import mean

# from xarray import corr


class ProcessIntrePro:
    """
    Process interpro data
    """

    def __init__(self, spark, data_file) -> None:
        self.spark = spark
        self.data_file = data_file
        self.spark_df = self.create_df(data_file)
        self.filter_df(self.spark_df)

    def create_df(self, file):
        """
        Create a pyspark data frame
        """
        spark_df = self.spark.read.csv(file, sep="\t", header=False, inferSchema=True)
        # col name info:
        # https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#output-formats
        col_names = [
            "protein_accession",
            "sequence_MD5",
            "sequence_length",
            "analysis",
            "signature_accession",
            "signature_description",
            "start_location",
            "stop_location",
            "score",
            "status",
            "date",
            "interpro_annot_accession",
            "interpro_annot_description",
            "GO",
            "pathway_annotations",
        ]
        spark_df = spark_df.toDF(*col_names)
        return spark_df

    def filter_df(self, spark_df):
        """
        Filter the data frame.
        """
        filtered_df = spark_df.filter(spark_df.interpro_annot_accession != "-")
        # filtered_df.select('interpro_annot_accession',).distinct().count()
        new_df = filtered_df.withColumn(
            "length_feature", filtered_df.stop_location - filtered_df.start_location
        )
        self.spark_df = new_df

    def distinct_protein_annots(self):
        """
        Get the distinct protein annotations.
        """
        n_distinct = self.spark_df.select("interpro_annot_accession").distinct()
        explain = n_distinct._sc._jvm.PythonSQLUtils.explainString(
            n_distinct._jdf.queryExecution(), "simple"
        )
        return n_distinct.count(), explain

    def get_number_annotations(self):
        """
        Get the number of annotations per protein accession.
        """
        n_annots = (
            self.spark_df.groupby("protein_accession")
            .agg(f.count("interpro_annot_accession"))
            .agg(f.mean("count(interpro_annot_accession)"))
        )
        explain = n_annots._sc._jvm.PythonSQLUtils.explainString(
            n_annots._jdf.queryExecution(), "simple"
        )
        n_annots = n_annots.collect()[0]["avg(count(interpro_annot_accession))"]

        return n_annots, explain

    def get_common_go_term(self):
        """
        Get the most common gene ontology (GO) terms.
        """
        gos_first = self.spark_df.select("GO").filter(
            (self.spark_df.GO != "-") & (~self.spark_df.GO.contains("|"))
        )

        gene_onts = self.spark_df.select(f.explode(f.split(self.spark_df.GO, "\|")))
        extra_gos = gene_onts.filter(gene_onts.col != "-")  # .show()

        all_go_terms = gos_first.union(extra_gos)

        common_gene_ont = all_go_terms.groupby("GO").count().orderBy(f.desc("count"))
        explain = common_gene_ont._sc._jvm.PythonSQLUtils.explainString(
            common_gene_ont._jdf.queryExecution(), "simple"
        )
        common_gene_ont = common_gene_ont.take(1)[0][0]

        return common_gene_ont, explain

    def get_avg_size_interpro(self):
        """
        Get the average size of an interpro.
        """
        avg_size = (
            self.spark_df.groupby("interpro_annot_accession")
            .agg(f.mean("length_feature"))
            .agg(f.mean("avg(length_feature)"))
        )
        explain = avg_size._sc._jvm.PythonSQLUtils.explainString(
            avg_size._jdf.queryExecution(), "simple"
        )
        avg_size = avg_size.collect()[0][0]

        return avg_size, explain

    def top_10_common_interpro(self):
        """
        Get the top 10 most common interpro accessions.
        """
        res_top_10 = (
            self.spark_df.groupby("interpro_annot_accession")
            .count()
            .orderBy(f.desc("count"))
        )
        explain = res_top_10._sc._jvm.PythonSQLUtils.explainString(
            res_top_10._jdf.queryExecution(), "simple"
        )
        res_top_10 = res_top_10.take(10)
        top_10 = [feature[0] for feature in res_top_10]

        return top_10, explain

    def top_10_big_features(self):
        """
        Get the top 10 biggest interpro features.
        """
        res_top_10_big = (
            self.spark_df.select("interpro_annot_accession")
            .filter(
                (
                    self.spark_df.length_feature
                    > (self.spark_df.sequence_length / 100 * 90)
                )
            )
            .groupby("interpro_annot_accession")
            .count()
            .orderBy(f.desc("count"))
        )
        explain = res_top_10_big._sc._jvm.PythonSQLUtils.explainString(
            res_top_10_big._jdf.queryExecution(), "simple"
        )
        res_top_10_big = res_top_10_big.take(10)
        top_10_big = [feature[0] for feature in res_top_10_big]

        return top_10_big, explain

    def top_10_words(self):
        """
        Get the top 10 most used words used in the interpro annotation
        description.
        """
        res_top_10_words = (
            self.spark_df.select(
                f.explode(
                    f.split(self.spark_df.interpro_annot_description, "\s+|,|-|\/")
                )
            )
            .filter(f.col("col") != "")
            .groupby(f.col("col"))
            .count()
            .orderBy(f.desc("count"))
        )
        explain = res_top_10_words._sc._jvm.PythonSQLUtils.explainString(
            res_top_10_words._jdf.queryExecution(), "simple"
        )
        res_top_10_words = res_top_10_words.take(10)
        top_10_words = [feature[0] for feature in res_top_10_words]

        return top_10_words, explain

    def top_10_least_common_words(self):
        """
        Get the top 10 least used words used in the interpro annotation
        description.
        """
        res_least_10_words = (
            self.spark_df.select(
                f.explode(
                    f.split(self.spark_df.interpro_annot_description, "\s+|,|-|\/")
                )
            )
            .filter(f.col("col") != "")
            .groupby(f.col("col"))
            .count()
            .orderBy(f.asc("count"))
        )
        explain = res_least_10_words._sc._jvm.PythonSQLUtils.explainString(
            res_least_10_words._jdf.queryExecution(), "simple"
        )
        res_least_10_words = res_least_10_words.take(10)
        least_10_words = [feature[0] for feature in res_least_10_words]

        return least_10_words, explain

    def common_words_large_features(self):
        """
        Get the most common words amongst large
        interpro features.
        """
        large_features = self.spark_df.filter(
            (self.spark_df.length_feature > (self.spark_df.sequence_length / 100 * 90))
        )
        res_top_10_words_large = (
            large_features.select(
                f.explode(
                    f.split(large_features.interpro_annot_description, "\s+|,|-|\/")
                )
            )
            .filter(f.col("col") != "")
            .groupby(f.col("col"))
            .count()
            .orderBy(f.desc("count"))
        )

        explain = res_top_10_words_large._sc._jvm.PythonSQLUtils.explainString(
            res_top_10_words_large._jdf.queryExecution(), "simple"
        )
        res_top_10_words_large = res_top_10_words_large.take(10)
        top_10_words_large = [feature[0] for feature in res_top_10_words_large]

        return top_10_words_large, explain

    def calc_corr(self):
        """
        Calculate the correlation.
        """
        correlation = self.spark_df.groupby("protein_accession").agg(
            f.mean("sequence_length"), f.count("interpro_annot_accession")
        )
        explain = correlation._sc._jvm.PythonSQLUtils.explainString(
            correlation._jdf.queryExecution(), "simple"
        )
        correlation = correlation.corr(
            "avg(sequence_length)", "count(interpro_annot_accession)"
        )

        return correlation, explain


def create_spark_sesson(n_threads):
    """
    Create a spark session
    """
    spark_context = SparkContext(f"local[{n_threads}]")
    spark = SparkSession(spark_context)
    return spark_context, spark


def make_data_dir(path: Path) -> None:
    """
    Create a directory (if it does not exsit yet) to store the
    data.

    :Excepts
    --------
    FileExistsError
        The directory already exists
    """
    try:
        path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print(f"[make_data_dir] {path} already exists.")


def main():
    """main"""
    output_dir = Path("output")
    make_data_dir(output_dir)
    n_threads = 16

    # Create a spark session
    spark_context, spark = create_spark_sesson(n_threads=n_threads)

    file = "/data/dataprocessing/interproscan/all_bacilli.tsv"

    process_interpro = ProcessIntrePro(spark, file)

    # Q1
    count, explain1 = process_interpro.distinct_protein_annots()
    # print(count)
    # print(explain1)
    n_annots, explain2 = process_interpro.get_number_annotations()
    # print(explain2)
    # print(n_annots)

    common_go, explain3 = process_interpro.get_common_go_term()
    # print(explain3)
    # print(common_GO)

    avg_size, explain4 = process_interpro.get_avg_size_interpro()
    # print(explain4)
    # print(avg_size)

    top_10, explain5 = process_interpro.top_10_common_interpro()
    # print(explain5)
    # print(top_10)

    top_10_big, explain6 = process_interpro.top_10_big_features()
    # print(explain6)
    # print(top_10_big)

    top_10_words, explain7 = process_interpro.top_10_words()
    # print(explain7)
    # print(top_10_words)

    top_10_least_words, explain8 = process_interpro.top_10_least_common_words()
    # print(explain8)
    # print(top_10_least_words)

    common_words_big_features, explain9 = process_interpro.common_words_large_features()
    # print(explain9)
    # print(common_words_big_features)

    correlation, explain10 = process_interpro.calc_corr()
    # print(explain10)
    # print(correlation)

    # field names
    fields = ["question", "answer", "explain"]

    questions = [i + 1 for i in range(10)]
    answers = [
        count,
        n_annots,
        common_go,
        avg_size,
        top_10,
        top_10_big,
        top_10_words,
        top_10_least_words,
        common_words_big_features,
        correlation,
    ]
    explanation = [
        explain1,
        explain2,
        explain3,
        explain4,
        explain5,
        explain6,
        explain7,
        explain8,
        explain9,
        explain10,
    ]

    rows = [[q, answ, expl] for q, answ, expl in zip(questions, answers, explanation)]

    # name of csv file
    filename = "assignment5.csv"

    # writing to csv file
    with open(output_dir / filename, "w", encoding="utf-8") as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile)

        # writing the fields
        csvwriter.writerow(fields)

        # writing the data rows
        csvwriter.writerows(rows)

    spark_context.stop()


if __name__ == "__main__":
    main()
