"""
Validation script for scOpAtlas structure

This script validates the module structure and code syntax without requiring
all dependencies to be installed.
"""

import os
import sys
import ast
from pathlib import Path

def check_file_syntax(filepath):
    """Check if a Python file has valid syntax"""
    try:
        with open(filepath, 'r') as f:
            code = f.read()
        ast.parse(code)
        return True, None
    except SyntaxError as e:
        return False, str(e)

def validate_module_structure():
    """Validate that all expected modules exist"""
    base_path = Path('/home/ubuntu/scidiff/scjdo/atlas')
    
    expected_files = [
        '__init__.py',
        'operator_metrics.py',
        'regime_classifier.py',
        'atlas_builder.py',
        'visualization.py',
        'build_atlas_cli.py'
    ]
    
    print("="*70)
    print("VALIDATING MODULE STRUCTURE")
    print("="*70)
    
    all_exist = True
    for filename in expected_files:
        filepath = base_path / filename
        exists = filepath.exists()
        status = "✅" if exists else "❌"
        print(f"{status} {filename}")
        
        if not exists:
            all_exist = False
    
    return all_exist

def validate_syntax():
    """Validate syntax of all Python files"""
    base_path = Path('/home/ubuntu/scidiff/scjdo/atlas')
    
    print("\n" + "="*70)
    print("VALIDATING PYTHON SYNTAX")
    print("="*70)
    
    all_valid = True
    for filepath in base_path.glob('*.py'):
        valid, error = check_file_syntax(filepath)
        status = "✅" if valid else "❌"
        print(f"{status} {filepath.name}")
        
        if not valid:
            print(f"   Error: {error}")
            all_valid = False
    
    return all_valid

def validate_imports():
    """Check that imports are properly structured"""
    print("\n" + "="*70)
    print("VALIDATING IMPORTS")
    print("="*70)
    
    # Check __init__.py
    init_file = Path('/home/ubuntu/scidiff/scjdo/atlas/__init__.py')
    
    try:
        with open(init_file, 'r') as f:
            code = f.read()
        
        tree = ast.parse(code)
        
        imports = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                if node.module:
                    for alias in node.names:
                        imports.append(alias.name)
        
        expected_imports = [
            'OperatorMetrics',
            'OperatorRegimeClassifier',
            'StableOperatorAtlas'
        ]
        
        all_present = True
        for imp in expected_imports:
            present = imp in imports
            status = "✅" if present else "❌"
            print(f"{status} {imp}")
            if not present:
                all_present = False
        
        return all_present
    
    except Exception as e:
        print(f"❌ Error checking imports: {e}")
        return False

def validate_classes():
    """Check that expected classes are defined"""
    print("\n" + "="*70)
    print("VALIDATING CLASS DEFINITIONS")
    print("="*70)
    
    files_and_classes = {
        'operator_metrics.py': ['OperatorMetrics'],
        'regime_classifier.py': ['OperatorRegimeClassifier'],
        'atlas_builder.py': ['StableOperatorAtlas']
    }
    
    base_path = Path('/home/ubuntu/scidiff/scjdo/atlas')
    all_valid = True
    
    for filename, expected_classes in files_and_classes.items():
        filepath = base_path / filename
        
        try:
            with open(filepath, 'r') as f:
                code = f.read()
            
            tree = ast.parse(code)
            
            defined_classes = []
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    defined_classes.append(node.name)
            
            for cls in expected_classes:
                present = cls in defined_classes
                status = "✅" if present else "❌"
                print(f"{status} {filename}: {cls}")
                if not present:
                    all_valid = False
        
        except Exception as e:
            print(f"❌ Error checking {filename}: {e}")
            all_valid = False
    
    return all_valid

def validate_methods():
    """Check that key methods are defined"""
    print("\n" + "="*70)
    print("VALIDATING KEY METHODS")
    print("="*70)
    
    # Check OperatorMetrics methods
    metrics_file = Path('/home/ubuntu/scidiff/scjdo/atlas/operator_metrics.py')
    
    expected_methods = {
        'OperatorMetrics': [
            'compute_jacobian',
            'compute_eigenvalues',
            'max_unstable_eigenvalue',
            'stability_depth',
            'plasticity_index',
            'stable_subspace_dim',
            'compute_all_metrics'
        ],
        'OperatorRegimeClassifier': [
            'classify',
            'get_regime_statistics',
            'compare_regimes_across_conditions'
        ],
        'StableOperatorAtlas': [
            'build',
            'validate_nonredundancy',
            'compare_conditions',
            'save'
        ]
    }
    
    all_valid = True
    
    for filename, class_methods in [
        ('operator_metrics.py', 'OperatorMetrics'),
        ('regime_classifier.py', 'OperatorRegimeClassifier'),
        ('atlas_builder.py', 'StableOperatorAtlas')
    ]:
        filepath = Path('/home/ubuntu/scidiff/scjdo/atlas') / filename
        
        try:
            with open(filepath, 'r') as f:
                code = f.read()
            
            tree = ast.parse(code)
            
            # Find the class
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef) and node.name == class_methods:
                    # Get method names
                    methods = [n.name for n in node.body if isinstance(n, ast.FunctionDef)]
                    
                    # Check expected methods
                    for method in expected_methods[class_methods]:
                        present = method in methods
                        status = "✅" if present else "❌"
                        print(f"{status} {class_methods}.{method}")
                        if not present:
                            all_valid = False
                    break
        
        except Exception as e:
            print(f"❌ Error checking {filename}: {e}")
            all_valid = False
    
    return all_valid

def validate_documentation():
    """Check that documentation files exist"""
    print("\n" + "="*70)
    print("VALIDATING DOCUMENTATION")
    print("="*70)
    
    docs = [
        'docs/scopatlas/DESIGN.md',
        'docs/scopatlas/README.md',
        'examples/tutorial_scopatlas.py'
    ]
    
    base_path = Path('/home/ubuntu/scidiff')
    all_exist = True
    
    for doc in docs:
        filepath = base_path / doc
        exists = filepath.exists()
        status = "✅" if exists else "❌"
        print(f"{status} {doc}")
        if not exists:
            all_exist = False
    
    return all_exist

def main():
    """Run all validations"""
    print("\n" + "="*70)
    print("scOpAtlas VALIDATION SUITE")
    print("="*70 + "\n")
    
    results = {
        'Module Structure': validate_module_structure(),
        'Python Syntax': validate_syntax(),
        'Imports': validate_imports(),
        'Class Definitions': validate_classes(),
        'Key Methods': validate_methods(),
        'Documentation': validate_documentation()
    }
    
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    all_passed = True
    for test, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status:10s} {test}")
        if not passed:
            all_passed = False
    
    print("="*70)
    
    if all_passed:
        print("\n🎉 All validations passed! scOpAtlas is ready to use.")
        return 0
    else:
        print("\n⚠️  Some validations failed. Please review the errors above.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
