#!/usr/bin/env python3
"""
Test cases for filter_probe.py
"""

import unittest
import tempfile
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.filter_probe import generate_kmers, process_probes, read_genome_kmers, get_cache_path, get_kmer_info

class TestFilterProbe(unittest.TestCase):
    def setUp(self):
        # Create temporary files for testing
        self.temp_files = []
        
    def tearDown(self):
        # Clean up temporary files
        for file in self.temp_files:
            if os.path.exists(file):
                os.remove(file)
                
    def create_temp_file(self, content):
        """Helper function to create temporary files with given content"""
        fd, path = tempfile.mkstemp()
        self.temp_files.append(path)
        with os.fdopen(fd, 'w') as tmp:
            tmp.write(content)
        return path

    def test_process_probes_basic(self):
        """Test basic probe processing functionality"""
        # Create a probe file with some sequences
        probe_content = """AAATTTCCC\tprobe1
GGGAAATTT\tprobe2
CCCGGGTTT\tprobe3"""
        probe_file = self.create_temp_file(probe_content)
        
        # Create genome k-mers (k=3)
        kmer_counts = {"AAA": 5, "TTT": 3, "CCC": 2}
        k = 3
        
        # Test processing
        processed = process_probes(probe_file, kmer_counts, k)
        self.assertEqual(len(processed), 3)
        # First probe should have AAA as max k-mer with count 5
        self.assertTrue("AAA\t5" in processed[0])
        
    def test_process_probes_empty_file(self):
        """Test with empty probe file"""
        probe_file = self.create_temp_file("")
        kmer_counts = {"AAA": 5, "TTT": 3}
        k = 3
        
        processed = process_probes(probe_file, kmer_counts, k)
        self.assertEqual(len(processed), 0)
        
    def test_process_probes_with_empty_lines(self):
        """Test handling of empty lines in probe file"""
        probe_content = """
AAATTTCCC\tprobe1

GGGAAATTT\tprobe2
"""
        probe_file = self.create_temp_file(probe_content)
        kmer_counts = {"AAA": 5}
        k = 3
        
        processed = process_probes(probe_file, kmer_counts, k)
        self.assertEqual(len(processed), 2)
        
    def test_process_probes_case_insensitive(self):
        """Test case-insensitive matching"""
        probe_content = """aaatttccc\tprobe1
AAATTTCCC\tprobe2"""
        probe_file = self.create_temp_file(probe_content)
        kmer_counts = {"AAA": 5}
        k = 3
        
        processed = process_probes(probe_file, kmer_counts, k)
        self.assertEqual(len(processed), 2)
        for line in processed:
            self.assertTrue("AAA\t5" in line)
            
    def test_get_kmer_info(self):
        """Test get_kmer_info function"""
        sequence = "AAATTTCCC"
        kmer_counts = {"AAA": 5, "TTT": 3, "CCC": 2}
        k = 3
        
        max_kmer, count = get_kmer_info(sequence, kmer_counts, k)
        self.assertEqual(max_kmer, "AAA")
        self.assertEqual(count, 5)
        
    def test_kmer_caching(self):
        """Test k-mer caching functionality"""
        # Create a test genome file
        genome_content = ">test\nAAAAATTTTTCCCCCGGGGG"
        genome_file = self.create_temp_file(genome_content)
        
        # First run - should compute and cache
        kmers1 = read_genome_kmers(genome_file, k=5)
        
        # Second run - should load from cache
        kmers2 = read_genome_kmers(genome_file, k=5)
        
        # Results should be identical
        self.assertEqual(kmers1, kmers2)
        
        # Check if cache file exists
        cache_file = get_cache_path(genome_file, k=5)
        self.assertTrue(cache_file.exists())
        
        # Add cache file to cleanup
        self.temp_files.append(str(cache_file))

if __name__ == '__main__':
    unittest.main()
