-- This file includes a couple of minor hacks which Django itself does
-- not manage. The SQL statements below should be run against the
-- repository database (in postgresql) once the initial "python
-- manage.py syncdb" step has completed.

--
-- This is a bit of a hack, the "library_extra" view is currently used
-- simply to convert library codes into their text plus numeric
-- components for the purposes of sorting in the web interface.
-- 
begin;

-- Simple-enough string to int casting, with handling for exceptions.
CREATE OR REPLACE FUNCTION convert_to_integer(v_input text)
RETURNS INTEGER AS $$
DECLARE v_int_value INTEGER DEFAULT NULL;
BEGIN
    BEGIN
        v_int_value := v_input::INTEGER;
    EXCEPTION WHEN OTHERS THEN
        RAISE NOTICE 'Invalid integer value: "%".  Returning NULL.', v_input;
        RETURN NULL;
    END;
RETURN v_int_value;
END;
$$ LANGUAGE plpgsql;

-- This is referenced in models.py as a
-- non-django-maintained table. It will be necessary to manage this
-- ourselves in e.g. testing fixtures.
create or replace view library_extra as
  select l.id as id,
  l.id as library_id,
  regexp_replace(l.code, '[0-9]+$', '') as code_text_prefix,
  convert_to_integer(regexp_replace(l.code, '^[a-zA-Z]+', '')) as code_numeric_suffix
  from library l;

-- Surprisingly, django doesn't seem to want to create constraints
-- enforcing table inheritance. We create them here. Note that this
-- block will fail and rollback on a schema which already possesses
-- these constraints (FIXME?).
begin;
alter table peakcalls add constraint peakcalls_dataprocess_fkey
  foreign key (dataprocess_ptr_id) references data_process(id) on delete cascade
  deferrable initially deferred;
alter table alignment add constraint alignment_dataprocess_fkey
  foreign key (dataprocess_ptr_id) references data_process(id) on delete cascade
  deferrable initially deferred;
alter table lane_qc add constraint lane_qc_dataprocess_fkey
  foreign key (dataprocess_ptr_id) references data_process(id) on delete cascade
  deferrable initially deferred;
commit;
