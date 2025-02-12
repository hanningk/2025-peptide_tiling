"""
===================
main
===================

Miscellaneous utilities for handling peptide tiles. v0.1
"""

#from typing import Union, Any, Optional

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from pandas import DataFrame

# Future:
    # pack the below into a class
    # inlude visualisation components

# Translate the input sequence.
# Modify this in future to better handle sequences indivisible by codon length
# and return that difference.

def initialise_sr(seqrecord: SeqRecord, offset_unit: int = 3, table: int | str = 1):
    """
    Validate and initialize a SeqRecord for translation.

    Args:
        seqrecord (SeqRecord): Input sequence record.
        offset_unit (int): Divisibility unit for sequence length (default: 3).
        table (int or str): Translation table (default: 1).

    Returns:
        tuple: The validated SeqRecord and its translated sequence (as Seq).

    Raises:
        ValueError: If the sequence is invalid or contains more than 1 stop codon.
        IndexError: If the sequence contains >1 stop codons or is not divisible by 3.
    """
    _check_offset_unit(seqrecord, offset_unit)

    tmp_trans_sr = seqrecord.seq.translate(table=table, to_stop=False)

    if len(seqrecord) % 3 == 0:

        if tmp_trans_sr.count("*") == 0:
            print(f"CDS: {seqrecord.seq}")
            print(f"Translation: {seqrecord.seq.translate(table=table, to_stop=True)}")
            return seqrecord, seqrecord.translate(table=table, to_stop=True)

        elif tmp_trans_sr.count("*") == 1:
            slice_end = len(seqrecord.translate(table=table, to_stop=True)) * 3
            seqrecord_trim = trim_sr(seqrecord, 0, slice_end)
            _check_offset_unit(seqrecord_trim, offset_unit)
            print(f"Trimming CDS to stop codon: {slice_end}")
            print(f"Trimmed CDS: {seqrecord_trim.seq}")
            print(f"Timmed Translation: {seqrecord_trim.seq.translate(table=table, to_stop=True)}")
            return seqrecord_trim, seqrecord_trim.translate(table=table, to_stop=True)

    else:
        raise IndexError(f"{seqrecord.id} is either indivisible by 3, or contains >1 stop codons.")

# Ensure sequence is divisible by offset unit
def _check_offset_unit(seq: SeqRecord | Seq, offset_unit: int) -> bool:
    """
    Validate that the sequence length is divisible by the offset unit.

    Args:
        seq (SeqRecord or Seq): The input sequence or SeqRecord.
        offset_unit (int): The unit to check divisibility.

    Returns:
        bool: True if divisible.

    Raises:
        ValueError: If the sequence length is not divisible by the offset unit.
    """
    length = len(seq.seq) if isinstance(seq, SeqRecord) else len(seq)
    if length % offset_unit == 0:
        return True
    raise ValueError(f"Sequence length {length} is not divisible by offset_unit {offset_unit}.")

def iterate_sr_to_tiles(seqrecord: SeqRecord,
    tile_len: int,
    tile_offset: int,
    offset_unit: int = 1) -> list[SeqRecord]:
    """
    Generates list of tiles as SeqRecords.

    Offset of final tile is reduced for consistent tile length.

    Args:
        seqrecord (SeqRecord): Sanitised SeqRecord.
        tile_len (int): Length of each tile.
        tile_offset (int): Offset of each subsequent tile.
        offset_unit (int): Handling of final tile offset, offset should be divisible by this.
    Returns:
        list[SeqRecord]:
    """
    # Initialize
    tile_start = 0-tile_offset
    tile_end = tile_len - tile_offset
    count = 0
    tiles = []

    # Generate SeqRecord for each tile
    while tile_end <= len(seqrecord):
        tile_start += tile_offset
        tile_end += tile_offset
        count += 1
        tile_record = trim_sr(seqrecord, tile_start, tile_end, count)

        if len(tile_record) == tile_len:
            tiles.append(tile_record)

        # Handle the final tile, modify the offset
        elif len(tile_record) < tile_len and tile_len - len(tile_record) < tile_offset:
            len_diff = tile_len - len(tile_record)
            final_offset = _final_tile_offset(len_diff, tile_offset, offset_unit)
            print(f"Offset of final tile reduced to {final_offset}")
            tile_start -= final_offset
            tile_end -= final_offset
            tile_record = trim_sr(seqrecord, tile_start, tile_end, count)
            tiles.append(tile_record)

            tile_start += tile_offset
            tile_end += tile_offset

    return tiles

def trim_sr(seqrecord: SeqRecord,
    slice_start: int,
    slice_end: int,
    count: int = None) -> SeqRecord:
    """Slice a SeqRecord into a tile, retaining features.

    Args:
        seqrecord (SeqRecord):
        slice_start (int):
        slice_end (int):
        count (int):
    Returns:
        (SeqRecord):
    """
    seqrecord_trim = seqrecord[slice_start:slice_end]
    if count:
        seqrecord_trim.name = f"{seqrecord.name}_partial_seq-tile_{count}"
        seqrecord_trim.id = f"{seqrecord.name}-slice_{count}"
        seqrecord_trim.description = f"{seqrecord.name}-tile_{count}-pos_{slice_start+1}-{slice_end}"
    else:
        seqrecord_trim.name = f"{seqrecord.name}_trimmed"
        seqrecord_trim.id = f"{seqrecord.name}-trimmed"
        seqrecord_trim.description = f"{seqrecord.name}-trimmed-pos_{slice_start+1}-{slice_end}"
    truncated_feats = get_overlap_feat(seqrecord, slice_start, slice_end)
    if truncated_feats:
        seqrecord_trim.features.extend(truncated_feats)
    return seqrecord_trim

def get_overlap_feat(seqrecord: SeqRecord,
    slice_start: int,
    slice_end: int) -> list[SeqFeature]:
    """Get features that paritally overlap tile bounds.

    Args:
        seqrecord (SeqRecord):
        slice_start (int):
        slice_end (int):

    Returns:
        (list[SeqFeature]):
    """
    truncated_feats = []
    for feat in seqrecord.features:

        feat_start = int(feat.location.start)
        feat_end = int(feat.location.end)

        if feat_start <= slice_end and slice_start <= feat_end and not (
            (slice_start <= feat_start and feat_end <= slice_end)):

            overlap_start = max(int(feat.location.start), slice_start)
            overlap_end = min(int(feat.location.end), slice_end)

            truncated_location = FeatureLocation(overlap_start, overlap_end)
            truncated_feat = feat.__class__(location=truncated_location, type=feat.type, qualifiers=feat.qualifiers.copy())

            truncated_feat.qualifiers.setdefault("note", []).append("truncated_feature")

            truncated_feats.append(truncated_feat)

    return truncated_feats

def _final_tile_offset(len_diff: int, tile_offset: int, offset_unit: int) -> int:
    if offset_unit == 0:
        raise ValueError("offset_unit cannot be zero")
    if type(offset_unit) != int:
        raise ValueError("offset_unit must be int type")
    if tile_offset % offset_unit == 0:
            return len_diff * offset_unit
    else:
        raise ValueError(f"tile_offset / offset_unit is not whole {tile_offset/offset_unit}")

## The following section is to convert tile reads to amino acid reads (dumb)

def tile_reads_to_aa(df: DataFrame,
    df_cols: list,
    tile_len: int,
    tile_offset: int,
    aa_info: DataFrame = None,
    aa_depth_col: str = 'depth',
    aa_feat_col: str = 'feat',
    aa_seq: str = None) -> DataFrame:

    # Initialise the output df
    num_rows = tile_len+(len(df)-1)*tile_offset
    df_out = DataFrame(0,index=range(num_rows), columns=df_cols, dtype='int64')

    # Check the provided depths
    if aa_info is None:
        print("Warning: aa_info is None, setting default depth to 1.")
        df_out[aa_depth_col] = 1
    elif len(df_out) != len(aa_info[aa_depth_col]):
        raise ValueError(f"Length mismatch: df_out has {len(df_out)} rows, "
                         f"but aa_info[{aa_depth_col}] has {len(aa_info[aa_depth_col])} rows. Exiting.")
    else:
        # Only executes if aa_info is not None and lengths match
        df_out[aa_depth_col] = aa_info[aa_depth_col]

    # Sum the tile counts across overlapping amino acids
    for i, row in df.iterrows():
        start_index = i * tile_offset
        end_index = start_index + tile_len
        if end_index > num_rows:
            end_index = num_rows
        values_to_add = row[df_cols].astype('int64').values
        df_out.loc[start_index:end_index - 1, df_cols] += values_to_add

    # Normalise per amino acid based on depth
    for k in df_cols:
        df_out[k+"_norm"] = df_out[k] / df_out[aa_depth_col]

    return df_out

def calc_percent(df: pd.DataFrame,
    df_cols: list):
    df_out = df
    columns_percent = []
    for i in df_cols:
        tot = df_out[i].sum()
        if tot != 0:
            df_out[i + '_percent'] = (df_out[i] / tot) * 100
        else:
            df_out[i + '_percent'] = 0
        columns_percent.append(i + '_percent')
    return df_out

def calc_enrichment(df: pd.DataFrame,
    pre_col: list,
    post_col: list):
    df_out = df
    for pre in pre_col:
        for post in post_col:
            enrich_col_name = f'enrichment_{post}_to_{pre}'

            # Calculate initial enrichment (this could give values >1 or <1)
            df_out[enrich_col_name] = df_out[post] / df_out[pre]
            # Adjust values: if enrichment is less than 1, take the inverse and make it negative
            df_out[enrich_col_name] = np.where(df_out[enrich_col_name] < 1,
                                                -(1 / df_out[enrich_col_name]),
                                                df[enrich_col_name])
    return df_out
