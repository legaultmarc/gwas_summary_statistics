-- Database creation script for PostgreSQL
BEGIN;


CREATE SCHEMA gwas_results;
COMMENT ON SCHEMA gwas_results IS 'Schema for the whole GWAS summary statistics management database.';


-- Types creation
CREATE TYPE gwas_results.statistical_test
    AS ENUM ('LogisticRegression', 'LinearRegression', 'CoxProportionalHazards',
             'SKAT', 'LMM_RepeatedMeasurements', 'LMM_PopulationStructure',
             'FisherExact');


-- Table creation
CREATE TABLE IF NOT EXISTS gwas_results.Variant (
    id serial,
    chrom varchar(2),
    position integer,
    reference varchar(100),
    alt varchar(100)[],
    CONSTRAINT pk_variant PRIMARY KEY (id)
);


CREATE TABLE IF NOT EXISTS gwas_results.GenotypingArray (
    name varchar(50),
    CONSTRAINT pk_geno_array PRIMARY KEY (name)
);


-- Many to many between Variant and GenotypingArray
CREATE TABLE IF NOT EXISTS gwas_results.ArrayVariants (
    variant integer REFERENCES gwas_results.Variant (id),
    genotyping_array varchar(50) REFERENCES gwas_results.GenotypingArray (name),
    CONSTRAINT pk_array_variants PRIMARY KEY (variant, genotyping_array)
);


CREATE TABLE IF NOT EXISTS gwas_results.Study (
    name varchar(50),
    parent varchar(50),
    genotyping_array varchar(50),
    CONSTRAINT pk_study PRIMARY KEY (name),
    CONSTRAINT fk_study_array FOREIGN KEY (genotyping_array) REFERENCES gwas_results.GenotypingArray (name),
    CONSTRAINT fk_study_study FOREIGN KEY (parent) REFERENCES gwas_results.Study (name)
);


CREATE TABLE IF NOT EXISTS gwas_results.StatModel (
    id serial,
    outcome varchar(50),
    subgroup varchar(50),
    covariates varchar(50)[],
    model gwas_results.statistical_test,
    software_name varchar(50),
    software_version varchar(20),
    n_samples integer,
    CONSTRAINT pk_statmodel PRIMARY KEY (id)
);


CREATE TABLE IF NOT EXISTS gwas_results.Results (
    model_id integer,
    p_value real,
    coefficient real,
    stderr real,
    variant_id integer,
    ref_allele varchar(100),
    imputed boolean,
    info real,
    mac integer,
    n integer,
    n_cases integer,
    n_controls integer,
    CONSTRAINT pk_gwas_results PRIMARY KEY (model_id, variant_id),
    CONSTRAINT fk_gr_model FOREIGN KEY (model_id) REFERENCES gwas_results.StatModel (id),
    CONSTRAINT fk_gr_variant FOREIGN KEY (variant_id) REFERENCES gwas_results.Variant (id)
);


-- Indexing
CREATE INDEX ON gwas_results.Variant(chrom, position, reference);


-- Functions
CREATE OR REPLACE FUNCTION gwas_results.InsertVariant(
    var_chrom varchar(2),
    var_pos integer,
    var_ref varchar(100),
    var_alt varchar(100),
    array_name varchar(50)
) RETURNS VOID AS
$$
DECLARE
    rec RECORD;
    last_variant integer;
BEGIN
    SELECT id, alt INTO STRICT rec FROM gwas_results.Variant WHERE chrom=var_chrom AND position=var_pos AND reference=var_ref;

    -- Is the alternative allele is in the array?
    IF NOT var_alt = ANY(rec.alt)
    THEN
        UPDATE gwas_results.Variant SET alt = rec.alt || var_alt WHERE id=rec.id;
    END IF;

    -- Is the variant in the many2many table for this array?
    IF NOT array_name = ANY(SELECT genotyping_array FROM gwas_results.ArrayVariants WHERE variant=rec.id)
    THEN
        INSERT INTO gwas_results.ArrayVariants (variant, genotyping_array)
            VALUES (rec.id, array_name);
    END IF;

    EXCEPTION
        WHEN NO_DATA_FOUND THEN
            -- The variant doesn't exist, so we add it
            INSERT INTO gwas_results.Variant (chrom, position, reference, alt[1])
                VALUES (var_chrom, var_pos, var_ref, var_alt)
                RETURNING id INTO last_variant;
            -- Inserting into the many2many table
            INSERT INTO gwas_results.ArrayVariants (variant, genotyping_array)
                VALUES (last_variant, array_name);

        WHEN TOO_MANY_ROWS THEN
            -- This is impossible...
            RAISE EXCEPTION 'Multiple entries for variant chr%:%/%', var_chrom, var_pos, var_ref;
END;
$$
LANGUAGE plpgsql STRICT;


CREATE OR REPLACE FUNCTION gwas_results.InsertVariantArray(
    var_chrom varchar(2)[],
    var_pos integer[],
    var_ref varchar(100)[],
    var_alt varchar(100)[],
    array_name varchar(50)
) RETURNS VOID AS
$$
BEGIN
    FOR i IN 1 .. array_upper(var_chrom, 1)
    LOOP
        PERFORM gwas_results.InsertVariant(var_chrom[i], var_pos[i], var_ref[i], var_alt[i], array_name);
    END LOOP;
END;
$$
LANGUAGE plpgsql STRICT;


COMMIT;
