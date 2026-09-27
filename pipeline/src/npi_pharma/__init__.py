"""NPI x pharmaceutical interaction scoring from an individual transcriptome.

Signature-composition engine (CMap / comboSC / DIPx style). Every score this
package emits is a Tier 3 (mechanism-only) plausibility ranking -- never a
measured or calibrated synergy.
"""

__version__ = "0.1.0"

EVIDENCE_TIER = "tier_3_mechanism_only"
