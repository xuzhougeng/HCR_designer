"""双探针设计配置模块"""
from dataclasses import dataclass, field
from pathlib import Path
import json
import logging

logger = logging.getLogger(__name__)

@dataclass
class DualProbeConfig:
    """双探针设计参数配置"""
    # 序列参数
    min_length: int = 15  
    max_length: int = 20
    gc_min: float = 40.0
    gc_max: float = 60.0
    tm_min: float = 45.0
    tm_max: float = 55.0
    
    # 间隔参数
    min_gap: int = 15  # 探针之间最小间距
    inner_gap: int = 0  # 探针内部间隔
    
    # 结构参数
    kmer_size: int = 8  # k-mer大小
    min_kmer_count: int = 2  # 最小k-mer计数
    min_complementary_length: int = 5  # 最小互补长度
    poly_n: int = 5  # 最大连续碱基数

    blast_db: str = None

    # 输出设置
    output_dir: Path = field(default=Path("output"))
    
    def __post_init__(self):
        """初始化后验证"""
        self.output_dir = Path(self.output_dir)
        self.validate()
    
    def validate(self):
        """验证参数有效性"""
        if self.min_length > self.max_length:
            raise ValueError(f"最小长度({self.min_length})不能大于最大长度({self.max_length})")
            
        if not 0 <= self.gc_min <= self.gc_max <= 100:
            raise ValueError(f"GC含量范围无效: {self.gc_min}-{self.gc_max}")
            
        if not 0 <= self.tm_min <= self.tm_max <= 100:
            raise ValueError(f"温度范围无效: {self.tm_min}-{self.tm_max}")
 
    def to_dict(self):
        """转换为字典格式"""
        return {
            "min_length": self.min_length,
            "max_length": self.max_length,
            "gc_min": self.gc_min,
            "gc_max": self.gc_max,
            "tm_min": self.tm_min,
            "tm_max": self.tm_max,
            "min_gap": self.min_gap,
            "inner_gap": self.inner_gap,
            "poly_n": self.poly_n,
            "kmer_size": self.kmer_size,
            "min_kmer_count": self.min_kmer_count,
            "output_dir": str(self.output_dir)
        }

    def save(self, filepath: Path):
        """保存配置到文件"""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)
