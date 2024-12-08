
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def analyze_blast_results(sequence, db_path):
    """
    对单个引物进行BLAST分析并统计不同错配数量的匹配数
    只考虑完全长度匹配的情况，gap也计入错配数

    Parameters:
    -----------
    sequence : str
        引物序列
    db_path : str
        BLAST数据库路径

    Returns:
    --------
    tuple: (dict, list)
        - dict: 包含不同错配数量的统计结果
        - list: 详细的匹配信息列表，每个元素为(subject_id, mismatches, gaps)
    """
    sequence_length = len(sequence)
    results = {}
    detailed_matches = []  # 存储详细的匹配信息
    
    # 创建临时文件存储查询序列
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp:
        temp.write(f">query\n{sequence}\n")
        query_file = temp.name

    try:
        # 运行BLAST
        output = NcbiblastnCommandline(
            query=query_file,
            db=db_path,
            outfmt="6 qseqid sseqid qstart qend sstart send nident length mismatch gapopen qseq sseq",
            word_size=10,
            task="blastn-short",
            dust="no",
            perc_identity=80
        )()[0]

        # 处理BLAST结果
        for line in output.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 12:
                continue

            # 解析BLAST结果字段
            subject_id = fields[1]     # 目标序列ID
            length = int(fields[7])    # 匹配长度
            mismatches = int(fields[8])# 错配数
            gaps = int(fields[9])      # gap数
            
            # 只考虑完全长度匹配的情况
            if length == sequence_length:
                # 计算总错配数（包括gaps）
                total_mismatches = mismatches + gaps
                
                # 更新统计结果
                if total_mismatches not in results:
                    results[total_mismatches] = 0
                results[total_mismatches] += 1
                
                # 保存详细匹配信息
                detailed_matches.append((subject_id, mismatches, gaps))

    finally:
        # 清理临时文件
        os.unlink(query_file)

    # 如果没有任何匹配结果，返回空结果
    if not results:
        return {0: 0}, []

    return results, detailed_matches

def get_detailed_blast_results(sequence, db_path):
    """
    获取详细的BLAST比对结果
    
    Parameters:
    -----------
    sequence : str
        引物序列
    db_path : str
        BLAST数据库路径
    
    Returns:
    --------
    list: 包含详细比对信息的记录列表
    """
    # 创建临时文件用于BLAST输入
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        f.write(f">query\n{sequence}\n")
        query_file = f.name

    # 设置BLAST参数
    blastn_cline = NcbiblastnCommandline(
        query=query_file,
        db=db_path,
        outfmt="6 qseqid sseqid qstart qend sstart send nident length mismatch gapopen qseq sseq",
        word_size=4,
        task="blastn-short",
        dust="no",
        perc_identity=80
    )

    # 运行BLAST
    stdout, stderr = blastn_cline()
    
    # 删除临时文件
    os.unlink(query_file)
    
    # 解析BLAST结果
    detailed_results = []
    
    for line in stdout.split('\n'):
        if not line.strip():
            continue
            
        fields = line.split('\t')
        if len(fields) < 12:
            continue
            
        subject_id = fields[1]
        mismatches = int(fields[8])
        query_seq = fields[10]
        subject_seq = fields[11]
        
        if mismatches <= 4:  # 只记录4个或更少错配的结果
            detailed_results.append({
                'subject_id': subject_id,
                'mismatches': mismatches,
                'alignment': f"Query:  {query_seq}\nSubject: {subject_seq}"
            })
    
    return detailed_results
