"""
Tests for FLYCOP.Consortia
Covers: initialization, species management, name lookup, common metabolite/reaction
        detection, and communications detection.
"""
import unittest
from unittest.mock import MagicMock

# FakeCobraModel is set up in conftest.py before this file is loaded
from conftest import FakeCobraModel
from FLYCOP.Consortia import Consortia


def _make_metabolite(met_id):
    m = MagicMock()
    m.id = met_id
    # exchange reaction for metabolites ending in _e
    rxn = MagicMock()
    rxn.id = "EX_" + met_id[:-2] if met_id.endswith("_e") else "other"
    m.reactions = [rxn]
    return m


def _make_model(model_id, metabolite_ids=None, reaction_ids=None):
    model = FakeCobraModel(model_id=model_id, name=model_id)
    model.metabolites = [_make_metabolite(m) for m in (metabolite_ids or [])]
    reaction_mocks = []
    for rid in (reaction_ids or []):
        rxn = MagicMock()
        rxn.id = rid
        reaction_mocks.append(rxn)
    model.reactions = reaction_mocks
    return model


class TestConsortiaInit(unittest.TestCase):
    def test_init_without_species_creates_empty_list(self):
        c = Consortia("test")
        self.assertEqual(c.name, "test")
        self.assertEqual(c.species, [])

    def test_init_with_valid_species(self):
        m = _make_model("ecoli")
        c = Consortia("test", species=[m])
        self.assertIn(m, c.species)

    def test_init_with_invalid_species_raises_type_error(self):
        with self.assertRaises(TypeError):
            Consortia("test", species=["not_a_model"])

    def test_class_level_species_does_not_leak_between_instances(self):
        """Mutable class-level list can cause leakage; init override prevents it."""
        c1 = Consortia("c1")
        c2 = Consortia("c2")
        c1.add_species(_make_model("ecoli"))
        self.assertEqual(len(c2.species), 0)


class TestConsortiaSpeciesManagement(unittest.TestCase):
    def setUp(self):
        self.c = Consortia("test")
        self.model_a = _make_model("A")
        self.model_b = _make_model("B")

    def test_add_species(self):
        self.c.add_species(self.model_a)
        self.assertIn(self.model_a, self.c.species)

    def test_remove_existing_species(self):
        self.c.add_species(self.model_a)
        self.c.remove_species(self.model_a)
        self.assertNotIn(self.model_a, self.c.species)

    def test_remove_nonexistent_species_does_not_raise(self):
        try:
            self.c.remove_species(self.model_a)
        except Exception as e:
            self.fail(f"remove_species raised unexpectedly: {e}")

    def test_get_species_returns_list(self):
        self.c.add_species(self.model_a)
        self.assertEqual(self.c.get_species(), [self.model_a])

    def test_get_species_names(self):
        self.c.add_species(self.model_a)
        self.c.add_species(self.model_b)
        names = self.c.get_species_names()
        self.assertIn("A", names)
        self.assertIn("B", names)

    def test_get_species_by_name_found(self):
        self.c.add_species(self.model_a)
        result = self.c.get_species_by_name("A")
        self.assertIs(result, self.model_a)

    def test_get_species_by_name_not_found_returns_none(self):
        result = self.c.get_species_by_name("nonexistent")
        self.assertIsNone(result)


class TestConsortiaUniformize(unittest.TestCase):
    def test_uniformize_accepts_self(self):
        """Bug was: missing self parameter — should not raise TypeError."""
        c = Consortia("test")
        try:
            c.uniformize()
        except TypeError as e:
            self.fail(f"uniformize() raised TypeError (likely missing self): {e}")


class TestConsortiaCommunications(unittest.TestCase):
    def test_get_common_metabolites(self):
        shared = _make_metabolite("glc_e")
        m1 = _make_model("A")
        m1.metabolites = [shared, _make_metabolite("atp_c")]
        m2 = _make_model("B")
        m2.metabolites = [shared, _make_metabolite("nadh_c")]

        c = Consortia("test", species=[m1, m2])
        common = c.get_common_metabolites("A", "B")
        self.assertIn(shared, common)

    def test_get_consortia_communications_detects_shared_exchange(self):
        shared_met = _make_metabolite("glc_e")  # ends with _e, has EX_ reaction
        m1 = _make_model("A")
        m1.metabolites = [shared_met]
        m2 = _make_model("B")
        m2.metabolites = [shared_met]

        c = Consortia("test", species=[m1, m2])
        comms = c.get_consortia_communications()
        self.assertEqual(len(comms), 1)
        self.assertEqual(comms[0], ("A", "B", "glc_e"))

    def test_get_consortia_communications_no_shared_returns_empty(self):
        m1 = _make_model("A")
        m1.metabolites = [_make_metabolite("atp_c")]  # no _e suffix
        m2 = _make_model("B")
        m2.metabolites = [_make_metabolite("nadh_c")]

        c = Consortia("test", species=[m1, m2])
        comms = c.get_consortia_communications()
        self.assertEqual(comms, [])


if __name__ == "__main__":
    unittest.main()
