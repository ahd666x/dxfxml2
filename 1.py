# این اسکریپت برای خواندن و پردازش فایل‌های DXF به کتابخانه 'ezdxf' نیاز دارد.
# لطفاً قبل از اجرا، این کتابخانه را نصب کنید (مثلاً با استفاده از: pip install ezdxf)

import ezdxf
import xml.etree.ElementTree as ET
import xml.dom.minidom # Import minidom for pretty printing
import math
import numpy as np
import re
import os
import sys
import subprocess # Import subprocess module
import unittest # Import the unittest framework

# --- Configuration ---
# Define DXF layer names used by the script
# This makes the script easier to update if layer names change
DXF_LAYER_CONFIG = {
    'part_border': '_ABF_PART_BORDER',
    'cutting_lines': '_ABF_CUTTING_LINES', # Still needed for front-only panels
    'sheet_border': '_ABF_SHEET_BORDER', # Layer for sheet borders
    'drilling_prefix': 'ABF_D', # Prefix for drilling layers (e.g., ABF_D10, ABF_D15)
    'pocket_layer': 'ABF_DSIDE_8', # Layer for pocket/side holes
    'groove_layer': 'ABF_GROOVE' # Layer for grooves
}

# --- Helper Functions ---

def get_bbox_center(vertices):
    """Calculates the center point of a bounding box from vertices."""
    if not vertices: return None, None
    min_x = min(v[0] for v in vertices)
    max_x = max(v[0] for v in vertices)
    min_y = min(v[1] for v in vertices)
    max_y = max(v[1] for v in vertices)
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    return center_x, center_y

def get_bbox(vertices):
    """Calculates the bounding box coordinates from vertices."""
    min_x = float('inf')
    max_x = -float('inf')
    min_y_val = float('inf')
    max_y_val = -float('inf')

    if not vertices: return min_x, min_y_val, max_x, max_y_val # Handle empty vertices

    for v in vertices:
        min_x = min(min_x, v[0])
        max_x = max(max_x, v[0])
        min_y_val = min(min_y_val, v[1])
        max_y_val = max(max_y_val, v[1])

    return min_x, min_y_val, max_x, max_y_val


def get_bbox_dimensions_sorted(vertices):
    """Calculates bounding box dimensions and returns them sorted (smallest, largest)."""
    min_x, min_y, max_x, max_y = get_bbox(vertices)
    if min_x == float('inf'): return 0.0, 0.0 # Handle empty vertices
    width = max_x - min_x
    height = max_y - min_y
    return round(min(width, height), 3), round(max(width, height), 3)

# Removed unused distance function
# def distance(p1, p2):
#     """Calculates the Euclidean distance between two 2D points."""
#     if p1 is None or p2 is None: return float('inf')
#     return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

def is_point_inside_bbox(point, bbox, tolerance=0.1):
    """Checks if a point is inside a bounding box with tolerance."""
    if point is None or bbox is None: return False
    return (bbox[0] - tolerance <= point[0] <= bbox[2] + tolerance and
            bbox[1] - tolerance <= point[1] <= bbox[3] + tolerance)


def get_number_from_layer_name_after_D(layer_name, drilling_prefix):
    """Extracts the number after the specified drilling prefix from layer names."""
    # Ensure the prefix is treated as literal text in the regex
    escaped_prefix = re.escape(drilling_prefix)
    match = re.search(rf'{escaped_prefix}(\d+)', layer_name.upper())
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return None
    return None

# --- Core Logic Functions ---

def find_and_group_panels(doc, config):
    """
    Finds _ABF_SHEET_BORDER, _ABF_PART_BORDER, and _ABF_CUTTING_LINES entities
    and groups them into physical panels based on containment within sheet borders
    and the presence of _ABF_PART_BORDER. Now only matches panels based on size,
    not position symmetry.
    Uses layer names from config.
    """
    sheet_border_layer = config['sheet_border'].upper()
    part_border_layer = config['part_border'].upper()
    cutting_lines_layer = config['cutting_lines'].upper()

    sheet_borders = [e for e in doc.modelspace() if e.dxftype() == 'LWPOLYLINE' and e.dxf.layer.upper() == sheet_border_layer]
    part_borders = [e for e in doc.modelspace() if e.dxftype() == 'LWPOLYLINE' and e.dxf.layer.upper() == part_border_layer]
    cutting_lines_borders = [e for e in doc.modelspace() if e.dxftype() == 'LWPOLYLINE' and e.dxf.layer.upper() == cutting_lines_layer]

    if len(sheet_borders) != 2:
        print(f"❌ Error: انتظار می شود دقیقاً دو موجودیت در لایه '{sheet_border_layer}' وجود داشته باشد، اما {len(sheet_borders)} یافت شد. نمی توان پنل ها را گروه بندی کرد.")
        return []

    # Determine which sheet border is "back" (contains part borders) and which is "front"
    sheet_border_bboxes = [get_bbox(list(sb.vertices())) for sb in sheet_borders]
    front_sheet_border_entity = None
    back_sheet_border_entity = None
    front_sheet_bbox = None
    back_sheet_bbox = None

    # Find the sheet border that contains part borders
    for i, sb_bbox in enumerate(sheet_border_bboxes):
        contains_part_border = False
        for pb in part_borders:
            pb_center = get_bbox_center(list(pb.vertices()))
            if pb_center and is_point_inside_bbox(pb_center, sb_bbox, tolerance=1.1):
                contains_part_border = True
                break
        if contains_part_border:
            back_sheet_border_entity = sheet_borders[i]
            back_sheet_bbox = sb_bbox
        else:
            front_sheet_border_entity = sheet_borders[i]
            front_sheet_bbox = sb_bbox

    if not front_sheet_border_entity or not back_sheet_border_entity:
         print(f"❌ Error: نتوانست مرزهای ورق جلو و پشت را بر اساس محتوای لایه '{part_border_layer}' شناسایی کند.")
         return []

    print(f"DEBUG: مرز ورق جلو شناسایی شد (محدوده: {front_sheet_bbox})")
    print(f"DEBUG: مرز ورق پشت شناسایی شد (محدوده: {back_sheet_bbox})")

    grouped_panels = []
    used_cutting_lines_indices = set()

    # 1. Process panels with _ABF_PART_BORDER (these are back-capable)
    for pb in part_borders:
        pb_vertices = list(pb.vertices())
        if not pb_vertices: continue
        pb_bbox = get_bbox(pb_vertices)
        pb_center = get_bbox_center(pb_vertices)

        # Confirm the part border is within the back sheet border
        if pb_center and is_point_inside_bbox(pb_center, back_sheet_bbox, tolerance=1.1):
            # Find matching cutting line border based ONLY on dimensions
            pb_dims = get_bbox_dimensions_sorted(pb_vertices)
            matched_cl_border = None
            matched_cl_index = -1

            # Search for a cutting line with matching dimensions
            for idx, cl_border in enumerate(cutting_lines_borders):
                 if idx not in used_cutting_lines_indices:
                      cl_vertices = list(cl_border.vertices())
                      if not cl_vertices: continue
                      cl_dims = get_bbox_dimensions_sorted(cl_vertices)
                      
                      # Check only for dimension match, ignore position
                      if math.isclose(pb_dims[0], cl_dims[0], abs_tol=0.1) and \
                         math.isclose(pb_dims[1], cl_dims[1], abs_tol=0.1):
                           matched_cl_border = cl_border
                           matched_cl_index = idx
                           print(f"DEBUG:   تطابق ابعاد پیدا شد برای مرز {part_border_layer} (ابعاد {pb_dims[0]:.1f} x {pb_dims[1]:.1f}) با مرز {cutting_lines_layer} (ابعاد {cl_dims[0]:.1f} x {cl_dims[1]:.1f}).")
                           break

            if matched_cl_border:
                 grouped_panels.append({
                     'borders': [pb, matched_cl_border],
                     'type': 'back_capable',
                     'primary_border': pb,
                     'secondary_border': matched_cl_border,
                     'primary_bbox': get_bbox(list(pb.vertices())),
                     'secondary_bbox': get_bbox(list(matched_cl_border.vertices())),
                     'sheet_border_front_bbox': front_sheet_bbox,
                     'sheet_border_back_bbox': back_sheet_bbox
                 })
                 used_cutting_lines_indices.add(matched_cl_index)
                 print(f"DEBUG: پنل back_capable گروه بندی شد با مرز ABF_CUTTING_LINES هم‌اندازه.")
            else:
                 print(f"WARNING: مرز ABF_PART_BORDER هم‌اندازه‌ای در لایه ABF_CUTTING_LINES یافت نشد.")
                 grouped_panels.append({
                     'borders': [pb],
                     'type': 'back_capable_standalone',
                     'primary_border': pb,
                     'secondary_border': None,
                     'primary_bbox': get_bbox(list(pb.vertices())),
                     'secondary_bbox': None,
                     'sheet_border_front_bbox': front_sheet_bbox,
                     'sheet_border_back_bbox': back_sheet_bbox
                 })
                 print(f"DEBUG: پنل back_capable_standalone گروه بندی شد.")

    # 2. Process remaining _ABF_CUTTING_LINES (front-only panels)
    for idx, cl_border in enumerate(cutting_lines_borders):
        if idx not in used_cutting_lines_indices:
            cl_vertices = list(cl_border.vertices())
            if not cl_vertices: continue
            cl_bbox = get_bbox(cl_vertices)
            cl_center = get_bbox_center(cl_vertices)

            if cl_center and is_point_inside_bbox(cl_center, front_sheet_bbox, tolerance=1.1):
                 grouped_panels.append({
                     'borders': [cl_border],
                     'type': 'front_only',
                     'primary_border': cl_border,
                     'secondary_border': None,
                     'primary_bbox': cl_bbox,
                     'secondary_bbox': None,
                     'sheet_border_front_bbox': front_sheet_bbox,
                     'sheet_border_back_bbox': back_sheet_bbox
                 })
                 print(f"DEBUG: مرز ABF_CUTTING_LINES به عنوان پنل فیزیکی تنها front_only اضافه شد.")
            else:
                 print(f"WARNING: مرز ABF_CUTTING_LINES خارج از محدوده ورق جلو یافت شد.")

    print(f"\nDEBUG: تعداد کل پنل‌های فیزیکی شناسایی شده پس از گروه بندی: {len(grouped_panels)}")
    return grouped_panels


def convert_coords_to_panel_system(point, panel_type, panel_xml_length, panel_xml_width,
                            primary_bbox_panel, secondary_bbox_panel,
                            sheet_border_front_bbox, sheet_border_back_bbox,
                            tolerance=0.1):
    """
    Converts DXF coordinates to XML panel coordinates with proper scaling and bounds checking.
    Returns a dictionary with transformed coordinates ('x', 'y') and determined face ('face').
    
    Args:
        point: Tuple of (x, y) coordinates in DXF space
        panel_type: Type of panel ('front_only' or 'back_capable')
        panel_xml_length: Length of the panel in XML coordinates
        panel_xml_width: Width of the panel in XML coordinates
        primary_bbox_panel: Tuple of (min_x, min_y, max_x, max_y) for primary panel bbox
        secondary_bbox_panel: Tuple of (min_x, min_y, max_x, max_y) for secondary panel bbox
        sheet_border_front_bbox: Tuple of (min_x, min_y, max_x, max_y) for front sheet bbox
        sheet_border_back_bbox: Tuple of (min_x, min_y, max_x, max_y) for back sheet bbox
        tolerance: Tolerance for coordinate checks (default: 0.1)
    
    Returns:
        Dictionary with 'x', 'y', and 'face' keys or None if conversion fails
    """
    if not isinstance(point, (list, tuple)) or len(point) < 2:
        print(f"DEBUG:   Invalid point format: {point}")
        return None
        
    entity_x, entity_y = point[0], point[1]
    print(f"DEBUG:   Converting DXF point ({entity_x:.3f}, {entity_y:.3f}) for panel (Type: {panel_type})")
    print(f"DEBUG:     Panel XML dimensions: {panel_xml_length:.0f} x {panel_xml_width:.0f}")

    # Determine which sheet the point is in
    in_front = is_point_inside_bbox((entity_x, entity_y), sheet_border_front_bbox, tolerance)
    in_back = is_point_inside_bbox((entity_x, entity_y), sheet_border_back_bbox, tolerance)
    
    # Select reference bbox and face based on sheet location
    reference_bbox = None
    face = None
    
    if in_front:
        reference_bbox = secondary_bbox_panel if panel_type == 'back_capable' else primary_bbox_panel
        face = "5"
    elif in_back:
        reference_bbox = primary_bbox_panel
        face = "6"
    else:
        print(f"DEBUG:     Point not in any sheet bbox")
        return None
        
    if not reference_bbox:
        print(f"DEBUG:     No valid reference bbox")
        return None

    # Get reference bbox dimensions and current point position
    ref_min_x, ref_min_y, ref_max_x, ref_max_y = reference_bbox
    ref_width = ref_max_x - ref_min_x
    ref_height = ref_max_y - ref_min_y
    
    # Calculate relative position within reference bbox (0-1)
    if ref_width == 0 or ref_height == 0:
        print(f"DEBUG:     Invalid reference bbox dimensions")
        return None
        
    rel_pos_x = (entity_x - ref_min_x) / ref_width
    rel_pos_y = (entity_y - ref_min_y) / ref_height
    
    # Scale to panel dimensions
    xml_x = rel_pos_x * panel_xml_length
    xml_y = rel_pos_y * panel_xml_width
    
    # Mirror Y coordinate for back face
    if face == "6":
        xml_y = panel_xml_width - xml_y
        
    # Ensure coordinates are within bounds
    xml_x = max(0.0, min(xml_x, panel_xml_length))
    xml_y = max(0.0, min(xml_y, panel_xml_width))
    
    print(f"DEBUG:     Transformed to XML coordinates: ({xml_x:.3f}, {xml_y:.3f}), Face: {face}")
    return {'x': xml_x, 'y': xml_y, 'face': face}


def create_drilling_xml(machines_element, entity, panel_type, panel_length, panel_width, 
                       primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                       sheet_border_back_bbox, tolerance, config):
    """Creates an XML element for a drilling operation."""
    if not hasattr(entity, 'dxf') or not hasattr(entity.dxf, 'center'):
        print(f"DEBUG:     خطا: موجودیت سوراخکاری فاقد مرکز است")
        return

    # Get the diameter from the layer name
    diameter = get_number_from_layer_name_after_D(entity.dxf.layer, config['drilling_prefix'])
    if diameter is None:
        print(f"DEBUG:     خطا: قطر سوراخکاری از نام لایه {entity.dxf.layer} قابل استخراج نیست")
        return

    # Extract center point
    center_point = (entity.dxf.center.x, entity.dxf.center.y)
    
    # Convert coordinates to panel system
    coords = convert_coords_to_panel_system(center_point, panel_type, panel_length, panel_width,
                                          primary_bbox_panel, secondary_bbox_panel,
                                          sheet_border_front_bbox, sheet_border_back_bbox)
    if coords is None:
        print(f"DEBUG:     خطا: تبدیل مختصات سوراخکاری ناموفق بود")
        return

    # Create the Machine element
    machine = ET.SubElement(machines_element, 'Machine')
    machine.set('Type', '2')  # Type 2 is for drilling
    machine.set('X', f"{coords['x']:.3f}")
    machine.set('Y', f"{coords['y']:.3f}")
    machine.set('Z', "0.000")  # Default Z position
    machine.set('Diameter', str(diameter))
    machine.set('Face', coords['face'])
    print(f"DEBUG:     سوراخکاری اضافه شد: مختصات=({coords['x']:.3f}, {coords['y']:.3f}), قطر={diameter}, سطح={coords['face']}")

def create_pocket_xml(machines_element, entity, panel_length, panel_width, panel_type,
                     primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                     sheet_border_back_bbox, tolerance, config):
    """Creates an XML element for a pocket operation."""
    if not hasattr(entity, 'get_points'):
        print(f"DEBUG:     خطا: موجودیت پاکت فاقد نقاط است")
        return

    vertices = [(p[0], p[1]) for p in entity.get_points()]
    if not vertices:
        print(f"DEBUG:     خطا: موجودیت پاکت فاقد نقطه است")
        return

    # Get center point of the pocket
    center_x, center_y = get_bbox_center(vertices)
    if center_x is None:
        print(f"DEBUG:     خطا: محاسبه مرکز پاکت ناموفق بود")
        return

    # Convert coordinates to panel system
    coords = convert_coords_to_panel_system((center_x, center_y), panel_type, panel_length, panel_width,
                                          primary_bbox_panel, secondary_bbox_panel,
                                          sheet_border_front_bbox, sheet_border_back_bbox)
    if coords is None:
        print(f"DEBUG:     خطا: تبدیل مختصات پاکت ناموفق بود")
        return

    # Calculate pocket dimensions
    min_x, min_y, max_x, max_y = get_bbox(vertices)
    width = abs(max_x - min_x)
    height = abs(max_y - min_y)

    # Create the Machine element for pocket
    machine = ET.SubElement(machines_element, 'Machine')
    machine.set('Type', '3')  # Type 3 is for pocket/side holes
    machine.set('X', f"{coords['x']:.3f}")
    machine.set('Y', f"{coords['y']:.3f}")
    machine.set('Z', "0.000")  # Default Z position
    machine.set('Width', f"{width:.3f}")
    machine.set('Height', f"{height:.3f}")
    machine.set('Face', coords['face'])
    print(f"DEBUG:     پاکت اضافه شد: مختصات=({coords['x']:.3f}, {coords['y']:.3f}), ابعاد=({width:.3f}x{height:.3f}), سطح={coords['face']}")

def create_groove_xml(machines_element, entity, panel_length, panel_width, panel_type,
                     primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                     sheet_border_back_bbox, tolerance, panel_thickness, config):
    """Creates an XML element for a groove operation."""
    if not hasattr(entity, 'get_points'):
        print(f"DEBUG:     خطا: موجودیت شیار فاقد نقاط است")
        return

    vertices = [(p[0], p[1]) for p in entity.get_points()]
    if not vertices:
        print(f"DEBUG:     خطا: موجودیت شیار فاقد نقطه است")
        return

    # Convert each vertex to panel coordinates
    converted_vertices = []
    for vertex in vertices:
        coords = convert_coords_to_panel_system(vertex, panel_type, panel_length, panel_width,
                                              primary_bbox_panel, secondary_bbox_panel,
                                              sheet_border_front_bbox, sheet_border_back_bbox)
        if coords is None:
            print(f"DEBUG:     خطا: تبدیل مختصات شیار ناموفق بود")
            return
        converted_vertices.append((coords['x'], coords['y']))
        face = coords['face']  # All vertices should have the same face

    # Create the Machine element for groove
    machine = ET.SubElement(machines_element, 'Machine')
    machine.set('Type', '4')  # Type 4 is for grooves
    
    # Add vertex points as subelements
    for i, (x, y) in enumerate(converted_vertices, 1):
        point = ET.SubElement(machine, f'Point{i}')
        point.set('X', f"{x:.3f}")
        point.set('Y', f"{y:.3f}")
        point.set('Z', "0.000")  # Default Z position

    machine.set('Face', face)
    machine.set('Depth', f"{panel_thickness/2:.3f}")  # Default groove depth is half panel thickness
    print(f"DEBUG:     شیار اضافه شد: تعداد نقاط={len(converted_vertices)}, سطح={face}")


def process_machining_entities_for_panel(doc, panel_element, panel_group_info, panel_length, panel_width, panel_thickness, config):
    """
    Iterates through all entities in the DXF model space, identifies
    those belonging to the current panel group, and processes them
    into Machining XML elements by calling appropriate helper functions.
    Uses layer names from config.
    """
    # Create or get the Machines element
    machines_element = panel_element.find('Machines')
    if machines_element is None:
        machines_element = ET.SubElement(panel_element, 'Machines')
    borders_in_group = panel_group_info['borders']
    panel_type = panel_group_info['type']
    # Get sheet border bboxes from panel_group_info
    sheet_border_front_bbox = panel_group_info['sheet_border_front_bbox']
    sheet_border_back_bbox = panel_group_info['sheet_border_back_bbox']
    # Get PANEL bboxes from panel_group_info
    primary_bbox_panel = panel_group_info['primary_bbox']
    secondary_bbox_panel = panel_group_info.get('secondary_bbox')

    # Calculate overall BBox for all borders in this panel group (for entity containment check)
    group_min_x, group_min_y = float('inf'), float('inf')
    group_max_x, group_max_y = -float('inf'), -float('inf')

    for border_entity_in_group in borders_in_group:
         border_vertices = list(border_entity_in_group.vertices())
         if not border_vertices: continue # Skip empty borders
         min_x, min_y, max_x, max_y = get_bbox(border_vertices)
         group_min_x = min(group_min_x, min_x)
         group_min_y = min(group_min_y, min_y)
         group_max_x = max(group_max_x, max_x)
         group_max_y = max(group_max_y, max_y)

    # Use the tolerance from convert_coords_to_panel_system (default 1.0) for containment check
    tolerance = 1.0 # Explicitly use the conversion tolerance here

    print(f"DEBUG: اسکان و پردازش موجودیت‌های ماشینکاری برای پنل فیزیکی...")

    for entity in doc.modelspace():
        # Debug message for every entity read from DXF
        if hasattr(entity, 'dxf') and hasattr(entity.dxf, 'layer'):
            print(f"DEBUG:   در حال بررسی موجودیت DXF: نوع={entity.dxftype()}, لایه={entity.dxf.layer}")
        else:
            print(f"DEBUG:   در حال بررسی موجودیت DXF: نوع={entity.dxftype()} (بدون اطلاعات لایه)")

        # Skip border entities themselves and sheet borders
        if entity in borders_in_group or entity.dxftype() == 'LWPOLYLINE' and entity.dxf.layer.upper() == config['sheet_border'].upper():
            print(f"DEBUG:     موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} یک مرز پنل یا مرز ورق است. رد می شود.")
            continue

        # Get the entity's reference point for checking containment within the panel group's overall BBox
        entity_point_to_check = None
        if hasattr(entity, 'dxf'):
             try:
                 if entity.dxftype() == 'CIRCLE':
                     entity_point_to_check = (entity.dxf.center.x, entity.dxf.center.y)
                 elif entity.dxftype() == 'LWPOLYLINE':
                     vertices = list(entity.vertices())
                     if vertices:
                         # Use the center of the polyline's bounding box
                         min_x, min_y, max_x, max_y = get_bbox(vertices)
                         entity_point_to_check = ((min_x + max_x) / 2.0, (min_y + max_y) / 2.0)
                 elif entity.dxftype() == 'LINE':
                     # Use midpoint of line
                     entity_point_to_check = ((entity.dxf.start.x + entity.dxf.end.x) / 2.0,
                                            (entity.dxf.start.y + entity.dxf.end.y) / 2.0)
                 elif entity.dxftype() in ['ARC', 'ELLIPSE']:
                     # Use center point
                     entity_point_to_check = (entity.dxf.center.x, entity.dxf.center.y)
                 elif entity.dxftype() == 'POINT':
                     entity_point_to_check = (entity.dxf.location.x, entity.dxf.location.y)
                 else:
                     print(f"DEBUG:     نوع موجودیت {entity.dxftype()} پشتیبانی نمی شود یا نقطه مرجع برای آن تعریف نشده است.")
             except Exception as e:
                 print(f"DEBUG:     خطا در استخراج نقطه مرجع برای موجودیت {entity.dxftype()}: {str(e)}")

        is_within_panel_group = False
        if entity_point_to_check is not None and len(entity_point_to_check) >= 2:
             # Check if the reference point is within either sheet border bbox
             is_within_front = is_point_inside_bbox(entity_point_to_check, sheet_border_front_bbox, tolerance)
             is_within_back = is_point_inside_bbox(entity_point_to_check, sheet_border_back_bbox, tolerance)
             is_within_panel_group = is_within_front or is_within_back
             print(f"DEBUG:     نقطه مرجع موجودیت: ({entity_point_to_check[0]:.3f}, {entity_point_to_check[1]:.3f})")
             print(f"DEBUG:     درون محدوده ورق جلو: {is_within_front}, درون محدوده ورق پشت: {is_within_back}")

        if not is_within_panel_group:
            print(f"DEBUG:     موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} خارج از محدوده پنل است. رد می شود.")
            continue

        # --- Process machining entities (only for entities within the overall BBox) ---

        # Process CIRCLEs (any layer starting with ABF_D) - Type="2".
        if entity.dxftype() == 'CIRCLE' and entity.dxf.layer.upper().startswith(config['drilling_prefix'].upper()):
            create_drilling_xml(machines_element, entity, panel_type, panel_length, panel_width, 
                              primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox, 
                              sheet_border_back_bbox, tolerance, config)
        elif entity.dxftype() == 'LWPOLYLINE' and entity.dxf.layer.upper() == config['pocket_layer'].upper():
            create_pocket_xml(machines_element, entity, panel_length, panel_width, panel_type,
                            primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                            sheet_border_back_bbox, tolerance, config)
        elif entity.dxftype() == 'LWPOLYLINE' and entity.dxf.layer.upper() == config['groove_layer'].upper():
            create_groove_xml(machines_element, entity, panel_length, panel_width, panel_type,
                            primary_bbox_panel, secondary_bbox_panel, sheet_border_front_bbox,
                            sheet_border_back_bbox, tolerance, panel_thickness, config)
        else:
            print(f"DEBUG:     موجودیت {entity.dxftype()} در لایه {entity.dxf.layer} با الگوهای ماشینکاری شناخته شده مطابقت ندارد. رد می شود.")


def create_panel_xml_structure(length, width, thickness):
    """Creates the basic XML structure for a panel"""
    # Create root element
    root = ET.Element('PanelAnalysis')
    
    # Add Panel element
    panel = ET.SubElement(root, 'Panel')
    panel.set('Length', f"{length:.3f}")
    panel.set('Width', f"{width:.3f}")
    panel.set('Thickness', f"{thickness:.3f}")
    
    # Add Machines element
    machines = ET.SubElement(panel, 'Machines')
    
    return root, panel


# --- Main Conversion Function ---

def dxf_to_custom_xml(input_file, config, panel_thickness=16.0): # Removed unused output_base_name
    """
    Main function to read DXF file, identify and process panels and their
    machining entities, and generate corresponding XML files.
    Uses layer names from config.
    """
    try:
        # Load the DXF document
        doc = ezdxf.readfile(input_file)
        print(f"DEBUG: فایل DXF '{input_file}' با موفقیت بارگذاری شد.")

        # 1. Find and group physical panels based on sheet borders and part borders
        grouped_panels = find_and_group_panels(doc, config)

        if not grouped_panels:
             print(f"❌ خطا: هیچ پنل فیزیکی برای پردازش یافت نشد.")
             return False # Indicate failure

        # Get the base name of the input DXF file without extension
        dxf_base_name = os.path.splitext(os.path.basename(input_file))[0]
        
        # Create output folder based on DXF file name
        output_folder = os.path.join(os.path.dirname(input_file), dxf_base_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            print(f"DEBUG: Created output folder: {output_folder}")

        # 2. Process each grouped physical panel
        for i, panel_group_info in enumerate(grouped_panels):
            primary_border_entity = panel_group_info['primary_border']
            # Determine XML panel dimensions (Length=larger, Width=smaller) from primary border
            # Note: get_bbox_dimensions_sorted returns (smallest, largest)
            # We need to assign the larger dimension to XML Length and smaller to XML Width
            panel_xml_width, panel_xml_length = get_bbox_dimensions_sorted(list(primary_border_entity.vertices()))
            length = panel_xml_length # Assign larger dimension to Length for XML
            width = panel_xml_width   # Assign smaller dimension to Width for XML

            print(f"\nDEBUG: --- شروع پردازش پنل فیزیکی شماره {i+1} (نوع: {panel_group_info['type']}) ---")
            print(f"DEBUG:   ابعاد پنل در XML (Length x Width): {length:.3f} x {width:.3f}")

            # 3. Create basic XML structure for the panel with the new naming convention
            # File name should be DXF_NAME.LengthxWidth.PanelIndex.xml for uniqueness
            output_file_name_base = f"{dxf_base_name}.{length:.0f}x{width:.0f}.{i+1}"
            output_file = os.path.join(output_folder, f"{output_file_name_base}.xml")

            # Panel ID and Name should be the same as the file name base (without .xml)
            panel_id = output_file_name_base
            panel_name = output_file_name_base


            root, panel_element = create_panel_xml_structure(panel_id, panel_name, length, width, panel_thickness)
            print(f"DEBUG:   ساختار اصلی XML برای پنل '{panel_name}' ایجاد شد.")

            # 4. Process machining entities within this panel group and add them to the XML structure
            # Pass all necessary info including panel bboxes and sheet bboxes
            process_machining_entities_for_panel(doc, panel_element, panel_group_info, length, width, panel_thickness, config)

            # 5. Save the generated XML file
            tree = ET.ElementTree(root)
            # Use indent for pretty printing the XML
            # Use minidom for pretty printing for compatibility with older Python versions
            xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
            dom = xml.dom.minidom.parseString(xml_string)
            # Ensure pretty_xml_string is bytes before writing to file in binary mode
            pretty_xml_string = dom.toprettyxml(indent="  ", encoding='utf-8')

            with open(output_file, "wb") as f:
                f.write(pretty_xml_string)


            print(f"✅ فایل '{output_file}' با موفقیت برای پنل فیزیکی شماره {i+1} (نوع: {panel_group_info['type']}) ایجاد شد.")
            print(f"DEBUG: --- پایان پردازش پنل فیزیکی شماره {i+1} ---")

        return True # Indicate success


    except FileNotFoundError:
        print(f"❌ خطا: فایل ورودی '{input_file}' یافت نشد.")
        return False # Indicate failure
    except ezdxf.DXFStructureError:
         print(f"❌ خطا: فایل '{input_file}' یک فایل DXF معتبر نیست یا خراب است.")
         return False # Indicate failure
    except Exception as e:
        print(f"❌ خطا در پردازش فایل DXF: {str(e)}")
        # Print full traceback for easier debugging of unexpected errors
        import traceback
        traceback.print_exc()
        return False # Indicate failure


def process_dxf_file(input_file, config):
    """
    Processes a DXF file to extract panel information and create XML files.
    Returns True if successful, False otherwise.
    """
    try:
        print(f"DEBUG: در حال خواندن فایل DXF: {input_file}")
        doc = ezdxf.readfile(input_file)
        
        # Find and group panels
        print("DEBUG: در حال جستجو و گروه‌بندی پنل‌ها...")
        panel_groups = find_and_group_panels(doc, config)
        if not panel_groups:
            print("❌ خطا: هیچ پنلی در فایل DXF یافت نشد.")
            return False

        print(f"DEBUG: {len(panel_groups)} پنل فیزیکی یافت شد.")

        # Process each panel group
        for i, panel_group_info in enumerate(panel_groups):
            print(f"\nDEBUG: --- شروع پردازش پنل فیزیکی شماره {i+1} ---")
            
            # Extract panel dimensions from the first border entity
            panel_dimensions = get_bbox_dimensions_sorted(panel_group_info['primary_border'].vertices())
            if not panel_dimensions or len(panel_dimensions) != 2:
                print(f"❌ خطا: استخراج ابعاد پنل شماره {i+1} ناموفق بود.")
                continue

            panel_width, panel_length = panel_dimensions
            panel_thickness = 18.0  # Default thickness in mm

            # Create output directory if it doesn't exist
            output_dir = os.path.splitext(input_file)[0]
            os.makedirs(output_dir, exist_ok=True)

            # Create output filename based on panel dimensions and number
            output_file = os.path.join(output_dir, 
                                     f"{os.path.splitext(os.path.basename(input_file))[0]}.{panel_length}x{panel_width}.{i+1}.xml")

            # Create XML structure
            print(f"DEBUG: ایجاد ساختار XML برای پنل {i+1}...")
            root, panel_element = create_panel_xml_structure(panel_length, panel_width, panel_thickness)
            
            # Get the Machines element for adding machining operations
            machines_element = panel_element.find('Machines')

            # Process machining entities for this panel
            print(f"DEBUG: پردازش موجودیت‌های ماشینکاری برای پنل {i+1}...")
            process_machining_entities_for_panel(doc, machines_element, panel_group_info,
                                               panel_length, panel_width, panel_thickness, config)

            # Save the XML file
            print(f"DEBUG: در حال ذخیره فایل XML برای پنل {i+1}...")
            # Use indent for pretty printing the XML
            xml_string = ET.tostring(root, encoding='utf-8').decode('utf-8')
            dom = xml.dom.minidom.parseString(xml_string)
            pretty_xml_string = dom.toprettyxml(indent="  ", encoding='utf-8')

            with open(output_file, "wb") as f:
                f.write(pretty_xml_string)

            print(f"✅ فایل '{output_file}' با موفقیت برای پنل فیزیکی شماره {i+1} (نوع: {panel_group_info['type']}) ایجاد شد.")
            print(f"DEBUG: --- پایان پردازش پنل فیزیکی شماره {i+1} ---")

        return True

    except FileNotFoundError:
        print(f"❌ خطا: فایل ورودی '{input_file}' یافت نشد.")
        return False
    except ezdxf.DXFStructureError:
        print(f"❌ خطا: فایل '{input_file}' یک فایل DXF معتبر نیست یا خراب است.")
        return False
    except Exception as e:
        print(f"❌ خطا در پردازش فایل DXF: {str(e)}")
        # Print full traceback for easier debugging of unexpected errors
        import traceback
        traceback.print_exc()
        return False


# --- Terminal UI Class ---
class TerminalUI:
    """Handles user interaction in the terminal."""
    def __init__(self, config):
        self.config = config

    def run(self):
        """Runs the terminal UI for selecting and converting a single DXF file."""
        while True: # Keep the UI running until explicitly exited
            try:
                subprocess.run(['clear' if os.name != 'nt' else 'cls'], shell=True)
            except:
                print('\n' * 100)
                
            print("=========================================")
            print("  مبدل DXF به XML برای پنل‌های چوبی")
            print("=========================================")

            dxf_files = self._get_dxf_files_in_current_directory()

            if not dxf_files:
                print("\n❌  هیچ فایل DXF در پوشه جاری یافت نشد.")
                print("لطفاً فایل‌های DXF را در کنار اسکریپت قرار دهید.")
                input("\nبرای خروج، Enter را فشار دهید...")
                sys.exit()

            print("\nفایل‌های DXF موجود در پوشه جاری:")
            for i, file_name in enumerate(dxf_files):
                print(f"{i + 1}. {file_name}")

            print("\n0. خروج")

            try:
                choice = int(input("\nلطفاً شماره فایل مورد نظر را وارد کنید (0-" + str(len(dxf_files)) + "): "))
                if 0 <= choice <= len(dxf_files):
                    if choice == 0:
                        print("\nخروج از برنامه.")
                        sys.exit()
                    else:
                        selected_file = dxf_files[choice - 1]
                        print(f"\n✅  فایل انتخابی: {selected_file}")
                        print("\nدر حال پردازش...")

                        # Run the conversion
                        success = dxf_to_custom_xml(selected_file, self.config)

                        if success:
                            print("\n✅  فرآیند تبدیل با موفقیت به پایان رسید.")
                        else:
                            print("\n❌  فرآیند تبدیل با مشکل مواجه شد.")

                        input("\nبرای ادامه، Enter را فشار دهید...")
                else:
                    print("⚠️  شماره نامعتبر.")
                    input("\nبرای ادامه، Enter را فشار دهید...")
            except ValueError:
                print("⚠️  ورودی نامعتبر.")
                input("\nبرای ادامه، Enter را فشار دهید...")

    def _get_dxf_files_in_current_directory(self):
        """Lists all files with .dxf extension in the current directory."""
        return [f for f in os.listdir('.') if os.path.isfile(f) and f.lower().endswith('.dxf')]


if __name__ == "__main__":
    print("DEBUG: شروع اجرای برنامه")
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        print(f"DEBUG: پردازش فایل ورودی: {input_file}")
        success = process_dxf_file(input_file, DXF_LAYER_CONFIG)
        sys.exit(0 if success else 1)
    else:
        # No command line argument provided, start the terminal UI
        ui = TerminalUI(DXF_LAYER_CONFIG)
        ui.run()
