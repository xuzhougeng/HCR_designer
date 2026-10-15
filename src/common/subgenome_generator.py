"""Build explicit subgenome annotations from FASTA headers and user-supplied features."""
from collections import Counter
import csv
import gzip
import io
from pathlib import Path
import re

from src.common.reference_annotation import MAX_LABEL_LENGTH


def default_label(filename):
    name = Path(filename).name
    if name.lower().endswith('.gz'):
        name = name[:-3]
    return re.sub(r'\s+', '_', Path(name).stem)


def validate_label(value):
    if not isinstance(value, str):
        raise ValueError('请填写亚基因组名称。')
    value = value.strip()
    if (not value or len(value) > MAX_LABEL_LENGTH or
            re.search(r'[\s,;"\x00-\x1f\x7f]', value)):
        raise ValueError(f'亚基因组名称需为 1–{MAX_LABEL_LENGTH} 个字符，不含空格、逗号、分号或双引号。')
    return value


def fasta_ids(upload):
    """Stream plain/gzip FASTA. Only the first word of each header is an identifier."""
    name = Path(upload.filename).name
    stem = name[:-3] if name.lower().endswith('.gz') else name
    if Path(stem).suffix.lower() not in ('.fa', '.fasta', '.fna'):
        raise ValueError(f'{name} 需为 .fa、.fasta、.fna 或其 .gz 压缩文件。')
    source = gzip.GzipFile(fileobj=upload.stream) if name.lower().endswith('.gz') else upload.stream
    found = False
    try:
        for number, raw in enumerate(source, 1):
            line = raw.decode('utf-8-sig').strip()
            if not line:
                continue
            if line.startswith('>'):
                parts = line[1:].split()
                if not parts:
                    raise ValueError(f'{name} 第 {number} 行缺少序列编号。')
                identifier = parts[0]
                # Keep CSV output readable by the existing two-column annotation loader.
                if (len(identifier) > 1024 or identifier.startswith('#') or
                        re.search(r'[,;"\x00-\x1f\x7f]', identifier)):
                    raise ValueError(f'{name} 第 {number} 行编号包含不支持的字符或超过 1024 位。')
                found = True
                yield identifier
            elif not found:
                raise ValueError(f'{name} 第 {number} 行的序列出现在 FASTA 标题之前。')
        if not found:
            raise ValueError(f'{name} 中没有 FASTA 序列编号。')
    except (OSError, EOFError, UnicodeError):
        raise ValueError(f'无法读取 {name}，请检查文件是否完整且为 UTF-8 FASTA。')
    finally:
        if source is not upload.stream:
            source.close()


def compile_rules(rules, position='contains', case_sensitive=False, rule_mode='feature'):
    if not isinstance(rules, list) or not 1 <= len(rules) <= 50:
        raise ValueError('请填写 1–50 条匹配规则。')
    if position not in ('contains', 'prefix', 'suffix') or rule_mode not in ('feature', 'literal'):
        raise ValueError('请选择有效的匹配方式和位置。')
    compiled = []
    for number, rule in enumerate(rules, 1):
        if not isinstance(rule, dict):
            raise ValueError(f'第 {number} 条规则格式无效。')
        pattern = rule.get('pattern')
        if not isinstance(pattern, str) or not pattern.strip() or len(pattern) > 200:
            raise ValueError(f'第 {number} 条规则需填写 1–200 个字符的编号特征。')
        pattern = pattern.strip()
        # Fixed-width digit placeholders; no arbitrary regular expressions or backtracking quantifiers.
        expression = ''.join('[0-9]' if rule_mode == 'feature' and char in 'xX' else re.escape(char)
                             for char in pattern)
        expression = ('^' if position == 'prefix' else '') + expression + ('$' if position == 'suffix' else '')
        compiled.append((re.compile(expression, 0 if case_sensitive else re.IGNORECASE),
                         validate_label(rule.get('label'))))
    return compiled


def csv_text(header, rows):
    output = io.StringIO(newline='')
    writer = csv.writer(output)
    writer.writerow(header)
    writer.writerows(rows)
    return output.getvalue()


def generate_subgenomes(uploads, mode, labels=None, rules=None, position='contains',
                       case_sensitive=False, rule_mode='feature'):
    if mode not in ('files', 'rules'):
        raise ValueError('请选择按文件或按编号规则生成。')
    if not uploads:
        raise ValueError('请先选择 FASTA 文件。')
    if mode == 'rules' and len(uploads) != 1:
        raise ValueError('按规则匹配时请选择一个 FASTA 文件。')
    if mode == 'files':
        if labels is None:
            labels = [default_label(upload.filename) for upload in uploads]
        if not isinstance(labels, list) or len(labels) != len(uploads):
            raise ValueError('每个 FASTA 文件都需要一个对应的亚基因组名称。')
        labels = [validate_label(label) for label in labels]
    compiled = compile_rules(rules, position, case_sensitive, rule_mode) if mode == 'rules' else []
    records, duplicate_count, header_count = {}, 0, 0
    for index, upload in enumerate(uploads):
        for identifier in fasta_ids(upload):
            header_count += 1
            folded = identifier.casefold()
            matches = {labels[index]} if mode == 'files' else {
                label for pattern, label in compiled if pattern.search(identifier)}
            if folded in records:
                duplicate_count += 1
                records[folded]['labels'].update(matches)
            else:
                records[folded] = dict(identifier=identifier, labels=matches)
    pairs, issues = [], []
    for record in records.values():
        identifier, matched = record['identifier'], sorted(record['labels'])
        if len(matched) == 1:
            pairs.append((identifier, matched[0]))
        else:
            issues.append(dict(identifier=identifier, issue='conflict' if matched else 'unmatched', labels=matched))
    conflicts = sum(issue['issue'] == 'conflict' for issue in issues)
    return dict(
        summary=dict(total=len(records), assigned=len(pairs), unmatched=len(issues) - conflicts,
                     conflicts=conflicts, duplicates=duplicate_count, headers=header_count,
                     files=len(uploads), subgenome_counts=dict(sorted(Counter(label for _, label in pairs).items()))),
        preview=[dict(identifier=identifier, subgenome=label) for identifier, label in pairs[:100]],
        issues=issues[:100], preview_limit=100,
        csv=csv_text(['transcript_id', 'subgenome'], pairs) if pairs and not conflicts else None,
        issues_csv=csv_text(['transcript_id', 'issue', 'subgenomes'],
                            [(issue['identifier'], issue['issue'], '|'.join(issue['labels'])) for issue in issues]) if issues else None)
