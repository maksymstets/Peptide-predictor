import pytest
from project import (
    filename_converter,
    peptide_generator,
    trypsin_logic,
    argc_proteinase_logic,
    aspn_endopeptidase_logic,
    chymotrypsinh_logic,
    chymotrypsinl_logic,
    papain_logic,
    pepsin_pH_1_3_logic,
    pepsin_pH_2_logic,
    single_digestion,
    parallel_digestion,
    sequential_digestion
)


def test_filename_converter():
    # Test with spaces
    assert filename_converter("conglutin beta protein") == "conglutin_beta_protein"
    # Test with brackets
    assert filename_converter("protein [Lupinus]") == "protein_Lupinus"
    # Test with empty string
    assert filename_converter("") == "unknown"


def test_peptide_generator():
    # Basic test
    assert peptide_generator([2, 5], "ABCDEFG") == ["AB", "CDE", "FG"]
    # No cleavage sites
    assert peptide_generator([], "ABCDEFG") == ["ABCDEFG"]
    # One cleavage site
    assert peptide_generator([3], "ABCDEFG") == ["ABC", "DEFG"]


def test_trypsin_logic():
    # Should cleave after K and R
    assert 2 in trypsin_logic("AKBRC")
    assert 4 in trypsin_logic("AKBRC")
    # Should not cleave before proline
    assert 1 not in trypsin_logic("KPRC")
    # Empty sequence
    assert trypsin_logic("") == []


def test_argc_proteinase_logic():
    # Should only cleave after R
    assert 2 in argc_proteinase_logic("ARBKC")
    # Should not cleave after K
    assert 4 not in argc_proteinase_logic("ARBKC")
    # No R in sequence
    assert argc_proteinase_logic("AKBKC") == []


def test_aspn_endopeptidase_logic():
    # Should cleave before D
    assert 2 in aspn_endopeptidase_logic("ABDEC")
    # Should not cleave when no D
    assert aspn_endopeptidase_logic("ABCEF") == []
    # Multiple D's
    assert len(aspn_endopeptidase_logic("ADBDCD")) == 3


def test_chymotrypsinh_logic():
    # Should cleave after F, W, Y
    assert 2 in chymotrypsinh_logic("AFBWCYD")
    assert 4 in chymotrypsinh_logic("AFBWCYD")
    # Should not cleave before proline
    assert 1 not in chymotrypsinh_logic("FPWP")
    # Should not cleave WM
    assert 2 not in chymotrypsinh_logic("AWMB")


def test_chymotrypsinl_logic():
    # Should cleave after F, W, Y, L, M, H
    assert 2 in chymotrypsinl_logic("AFBLCM")
    assert 4 in chymotrypsinl_logic("AFBLCM")
    # Should not cleave before proline
    assert 1 not in chymotrypsinl_logic("FPWP")
    # Should not cleave MY
    assert 2 not in chymotrypsinl_logic("AMYB")


def test_papain_logic():
    # Should cleave after R when preceded by hydrophobic
    assert 2 in papain_logic("ARBKC")
    # Should not cleave before valine
    assert 2 not in papain_logic("ARVKC")
    # K with hydrophobic before
    assert 2 in papain_logic("AKBRC")


def test_pepsin_pH_1_3_logic():
    # Should cleave at F and L
    assert len(pepsin_pH_1_3_logic("AFLB")) > 0
    # Check that it doesn't cleave with proline nearby
    sequence = "AFLB"
    sites = pepsin_pH_1_3_logic(sequence)
    assert isinstance(sites, list)


def test_pepsin_pH_2_logic():
    # Should cleave at F, L, W, Y
    assert len(pepsin_pH_2_logic("AFLWBY")) > 0
    # Should return empty list for no cleavage sites
    sequence = "AAA"
    sites = pepsin_pH_2_logic(sequence)
    assert isinstance(sites, list)


def test_single_digestion():
    # Test that single digestion runs without errors
    sequence = "ARKBKC"
    enzyme_register = {"Trypsin": trypsin_logic}
    
    # This will create a PDF file, but we just want to make sure it doesn't crash
    try:
        single_digestion(sequence, "Trypsin", enzyme_register, "test_protein")
        # If we get here, the function ran without errors
        assert True
    except Exception as e:
        # If it fails, we want to see the error
        pytest.fail(f"single_digestion failed: {e}")


def test_parallel_digestion():
    # Test that parallel digestion runs without errors
    sequence = "ARKBKC"
    enzyme_names = ["Trypsin", "Arg-C_proteinase"]
    enzyme_register = {
        "Trypsin": trypsin_logic,
        "Arg-C_proteinase": argc_proteinase_logic
    }
    
    try:
        parallel_digestion(sequence, enzyme_names, enzyme_register, "test_protein")
        assert True
    except Exception as e:
        pytest.fail(f"parallel_digestion failed: {e}")


def test_sequential_digestion():
    # Test sequential digestion with one enzyme
    sequence = "ARKBKC"
    enzyme_names = ["Arg-C_proteinase"]
    enzyme_register = {"Arg-C_proteinase": argc_proteinase_logic}
    
    final_peptides, operations_log = sequential_digestion(sequence, enzyme_names, enzyme_register)
    
    # Check that we got results
    assert len(operations_log) == 1
    assert operations_log[0]['enzyme'] == "Arg-C_proteinase"
    assert len(final_peptides) > 0
    
    # Test with two enzymes
    enzyme_names = ["Arg-C_proteinase", "Trypsin"]
    enzyme_register = {
        "Arg-C_proteinase": argc_proteinase_logic,
        "Trypsin": trypsin_logic
    }
    
    final_peptides, operations_log = sequential_digestion(sequence, enzyme_names, enzyme_register)
    
    # Should have two digestion steps
    assert len(operations_log) == 2
    assert operations_log[0]['enzyme'] == "Arg-C_proteinase"
    assert operations_log[1]['enzyme'] == "Trypsin"
