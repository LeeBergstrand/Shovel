#!/usr/bin/env python

from Trim_SPAdes_FASTA import FASTAFormat, filter_seq_object
from unittest import TestCase


class TestYAMLRead(TestCase):
    def test_coverage_filtering(self):
        """Tests that filtering function is properly filtering by coverage."""
        extra_low_coverage = FASTAFormat(20, 2000, 0.00, 80, 0, "ATCG")
        self.assertTrue(filter_seq_object(extra_low_coverage, 1000, 80.00))

        low_coverage = FASTAFormat(20, 2000, 79.00, 81, 0, "ATCG")
        self.assertTrue(filter_seq_object(low_coverage, 1000, 80.00))

        equal_coverage = FASTAFormat(20, 2000, 80.00, 82, 0, "ATCG")
        self.assertFalse(filter_seq_object(equal_coverage, 1000, 80.00))

        high_coverage = FASTAFormat(20, 2000, 81.00, 83, 0, "ATCG")
        self.assertFalse(filter_seq_object(high_coverage, 1000, 80.00))

        extra_high_coverage = FASTAFormat(20, 2000, 100.00, 84, 0, "ATCG")
        self.assertFalse(filter_seq_object(extra_high_coverage, 1000, 80.00))

    def test_length_filtering(self):
        """Tests that filtering function is properly filtering by length."""
        extra_low_coverage = FASTAFormat(20, 0, 81.00, 85, 0, "ATCG")
        self.assertTrue(filter_seq_object(extra_low_coverage, 1000, 80.00))

        low_coverage = FASTAFormat(20, 999, 81.00, 86, 0, "ATCG")
        self.assertTrue(filter_seq_object(low_coverage, 1000, 80.00))

        equal_coverage = FASTAFormat(20, 1000, 81.00, 87, 0, "ATCG")
        self.assertFalse(filter_seq_object(equal_coverage, 1000, 80.00))

        high_coverage = FASTAFormat(20, 1001, 81.00, 88, 0, "ATCG")
        self.assertFalse(filter_seq_object(high_coverage, 1000, 80.00))

        extra_high_coverage = FASTAFormat(20, 2000, 81.00, 89, 0, "ATCG")
        self.assertFalse(filter_seq_object(extra_high_coverage, 1000, 80.00))
