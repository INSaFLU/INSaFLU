"""
pytest configuration for BFS diagnostic tests.
"""

import os
import sys

# Add the INSaFLU directory to Python path for imports
INSAFLU_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, INSAFLU_DIR)

