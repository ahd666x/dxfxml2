import os
import xml.etree.ElementTree as ET

xml_dir = '1'
for xml_file in sorted(os.listdir(xml_dir)):
    if not xml_file.endswith('.xml'):
        continue
        
    full_path = os.path.join(xml_dir, xml_file)
    print(f"\n=== {xml_file} ===")
    try:
        with open(full_path, 'rb') as f:
            content = f.read()
            print(f'File size: {len(content)} bytes')
        
        tree = ET.parse(full_path)
        root = tree.getroot()
        
        # Print basic structure
        print("\nStructure:")
        def print_element(elem, level=0):
            print("  " * level + f"<{elem.tag} {' '.join(f'{k}=\"{v}\"' for k,v in elem.attrib.items())}>")
            for child in elem:
                print_element(child, level + 1)
                
        print_element(root)
    except Exception as e:
        print(f'Error: {e}')

