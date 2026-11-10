"""Validate the current fresh pilot and its reusable profiles, without a historic cache."""
import sys
sys.dont_write_bytecode = True
from verify_rebuild import main

if __name__ == "__main__":
    main()
