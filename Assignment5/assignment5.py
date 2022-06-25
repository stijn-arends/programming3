import pyspark
from pyspark import SparkContext
from pyspark import SparkFiles
from pyspark.sql import SQLContext
from pyspark.sql import SparkSession
from pyspark.sql.functions import mean
import pyspark.sql.functions as f
import csv
from pathlib import Path

# from xarray import corr



class ProcessIntrePro:

    def __init__(self, spark, data_file):
        self.spark = spark
        self.data_file = data_file
        self.df = self.create_df(data_file)
        self.filter_df(self.df)

    def create_df(self, file):
        df = self.spark.read.csv(file, sep="\t", header=False, inferSchema=True)
        # col name info: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#output-formats
        col_names = ["protein_accession", "sequence_MD5", "sequence_length", 
            "analysis", "signature_accession", "signature_description", "start_location",
            "stop_location", "score", "status", "date", "interpro_annot_accession", "interpro_annot_description", 
            "GO", "pathway_annotations"]
        df = df.toDF(*col_names)
        return df

    def filter_df(self, df):
        filtered_df = df.filter(df.interpro_annot_accession != "-")
        # filtered_df.select('interpro_annot_accession',).distinct().count()  
        new_df = filtered_df.withColumn('length_feature',
                       filtered_df.stop_location - filtered_df.start_location)
        self.df = new_df

    def distinct_protein_annots(self):
        n_distinct = self.df.select('interpro_annot_accession').distinct()
        explain = n_distinct._sc._jvm.PythonSQLUtils.explainString(n_distinct._jdf.queryExecution(), "simple")
        return n_distinct.count(), explain

    
    def get_number_annotations(self):
        n_annots = self.df.groupby('protein_accession').agg(f.count('interpro_annot_accession')).agg(f.mean('count(interpro_annot_accession)'))
        explain = n_annots._sc._jvm.PythonSQLUtils.explainString(n_annots._jdf.queryExecution(), "simple")
        n_annots = n_annots.collect()[0]['avg(count(interpro_annot_accession))']

        return n_annots, explain


    def get_common_GO_term(self):
        gos_first = self.df.select('GO').filter((self.df.GO != "-") & (~self.df.GO.contains('|')))

        GOs = self.df.select(f.explode(f.split(self.df.GO, '\|')))
        extra_gos = GOs.filter(GOs.col != "-") # .show()

        all_go_terms = gos_first.union(extra_gos)

        common_GO = all_go_terms.groupby('GO').count().orderBy(f.desc('count'))
        explain = common_GO._sc._jvm.PythonSQLUtils.explainString(common_GO._jdf.queryExecution(), "simple")
        common_GO = common_GO.take(1)[0][0]
        
        return common_GO, explain


    def get_avg_size_interpro(self):
        avg_size = self.df.groupby('interpro_annot_accession').agg(f.mean('length_feature')).agg(f.mean('avg(length_feature)'))
        explain = avg_size._sc._jvm.PythonSQLUtils.explainString(avg_size._jdf.queryExecution(), "simple")
        avg_size = avg_size.collect()[0][0]

        return avg_size, explain


    def top_10_common_interpro(self):
        res_top_10 = self.df.groupby('interpro_annot_accession').count().orderBy(f.desc('count'))
        explain = res_top_10._sc._jvm.PythonSQLUtils.explainString(res_top_10._jdf.queryExecution(), "simple")
        res_top_10 = res_top_10.take(10)
        top_10 = [feature[0] for feature in res_top_10]
        
        return top_10, explain

    def top_10_big_features(self):
        res_top_10_big = self.df.select('interpro_annot_accession').filter((self.df.length_feature > (self.df.sequence_length / 100 * 90))).groupby('interpro_annot_accession').count().orderBy(f.desc('count'))
        explain = res_top_10_big._sc._jvm.PythonSQLUtils.explainString(res_top_10_big._jdf.queryExecution(), "simple")
        res_top_10_big = res_top_10_big.take(10)
        top_10_big = [feature[0] for feature in res_top_10_big]
        
        return top_10_big, explain


    def top_10_words(self):
        res_top_10_words = self.df.select(f.explode(f.split(self.df.interpro_annot_description, '\s+|,|-|\/'))).filter(f.col("col") != "").groupby(f.col('col')).count().orderBy(f.desc('count'))
        explain = res_top_10_words._sc._jvm.PythonSQLUtils.explainString(res_top_10_words._jdf.queryExecution(), "simple")
        res_top_10_words = res_top_10_words.take(10)
        top_10_words = [feature[0] for feature in res_top_10_words]
        
        return top_10_words, explain

    def top_10_least_common_words(self):
        res_least_10_words = self.df.select(f.explode(f.split(self.df.interpro_annot_description, '\s+|,|-|\/'))).filter(f.col("col") != "").groupby(f.col('col')).count().orderBy(f.asc('count'))
        explain = res_least_10_words._sc._jvm.PythonSQLUtils.explainString(res_least_10_words._jdf.queryExecution(), "simple")
        res_least_10_words = res_least_10_words.take(10)
        least_10_words = [feature[0] for feature in res_least_10_words]
        
        return least_10_words, explain

    def common_words_large_features(self):
        large_features = self.df.filter((self.df.length_feature > (self.df.sequence_length / 100 * 90)))
        res_top_10_words_large = large_features.select(f.explode(f.split(large_features.interpro_annot_description, '\s+|,|-|\/'))).filter(f.col("col") != "").groupby(f.col('col')).count().orderBy(f.desc('count'))

        explain = res_top_10_words_large._sc._jvm.PythonSQLUtils.explainString(res_top_10_words_large._jdf.queryExecution(), "simple")
        res_top_10_words_large = res_top_10_words_large.take(10)
        top_10_words_large = [feature[0] for feature in res_top_10_words_large]
       
        return top_10_words_large, explain

    def calc_corr(self):
        correlation = self.df.groupby('protein_accession').agg(f.mean('sequence_length'), f.count('interpro_annot_accession'))
        explain = correlation._sc._jvm.PythonSQLUtils.explainString(correlation._jdf.queryExecution(), "simple")
        correlation = correlation.corr("avg(sequence_length)", 'count(interpro_annot_accession)')

        return correlation, explain




def create_spark_sesson(n_threads):
    sc = SparkContext(f'local[{n_threads}]')
    spark = SparkSession(sc)
    return sc, spark


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

    output_dir = Path("output")
    make_data_dir(output_dir)
    n_threads = 16

    # Create a spark session
    sc, spark = create_spark_sesson(n_threads=n_threads)

    file = "/data/dataprocessing/interproscan/all_bacilli.tsv"

    process_interpro = ProcessIntrePro(spark, file)

    # Q1
    count, explain1 = process_interpro.distinct_protein_annots()
    # print(count)
    # print(explain1)
    n_annots, explain2 = process_interpro.get_number_annotations()
    # print(explain2)
    # print(n_annots)

    common_GO, explain3 = process_interpro.get_common_GO_term()
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
    fields = ['question', 'answer', 'explain']

    questions = [i+1 for i in range(10)] 
    answers = [count, n_annots, common_GO, avg_size, top_10, top_10_big, top_10_words, 
        top_10_least_words, common_words_big_features, correlation]
    explanation = [explain1, explain2, explain3, explain4, explain5, explain6, explain7, explain8,
        explain9, explain10]

    rows = [[q, answ, expl] for q, answ, expl in zip(questions, answers, explanation)]
        
    # name of csv file 
    filename = "assignment5.csv"
        
    # writing to csv file 
    with open(output_dir / filename, 'w') as csvfile: 
        # creating a csv writer object 
        csvwriter = csv.writer(csvfile) 
            
        # writing the fields 
        csvwriter.writerow(fields) 
            
        # writing the data rows 
        csvwriter.writerows(rows)

    sc.stop()


if __name__ == "__main__":
    main()