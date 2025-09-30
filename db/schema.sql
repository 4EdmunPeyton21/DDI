CREATE TABLE drugs (
  drugbank_id TEXT PRIMARY KEY,
  rxnorm_id TEXT,
  name TEXT,
  synonyms TEXT[],
  atc_codes TEXT[],
  smiles TEXT,
  properties JSONB,
  raw JSONB,
  created_at TIMESTAMP DEFAULT now()
);

CREATE TABLE interactions (
  id SERIAL PRIMARY KEY,
  drug_a TEXT NOT NULL,
  drug_b TEXT NOT NULL,
  severity TEXT,
  mechanism TEXT,
  effects TEXT[],
  management TEXT,
  evidence JSONB,
  source TEXT,
  updated_at TIMESTAMP DEFAULT now(),
  UNIQUE (drug_a, drug_b)
);